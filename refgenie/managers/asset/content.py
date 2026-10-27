"""
AssetContentManager - the content write path, plus its stateless helpers.

Reached as ``rgc.asset.content``. It turns produced content (a build output, a
pulled archive, a folder the user registers) into catalog rows: content
addressing, placement under ``data/``, and the name rows. Its reads go straight
to ``queries``, so it holds no reference to ``AssetManager``.

The module functions run inside a caller's open session or use no manager
state at all.
"""

import os
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any, TYPE_CHECKING

from sqlalchemy.engine import Engine
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.tables import Asset, AssetGroup, AssetName, Genome, Recipe, SeekKeyType
from refgenie.exceptions import MissingAssetError, MissingAssetGroupError, MissingGenomeError
from refgenie.logger import logger
from refgenie.managers.asset.group import asset_group_in_session
from refgenie.managers.asset.queries import asset_by_digest_stmt, asset_by_name_stmt
from refgenie.managers.asset.seek_key import bind_asset_class_seek_keys
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.models import GenomeDigest
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import Change
from refgenie.utils.digest import BuildProvenance, directory_size, get_dir_digest

if TYPE_CHECKING:
    from refgenie.managers.asset.alias_tree import AliasTree
    from refgenie.managers.asset.group import AssetGroupManager
    from refgenie.managers.asset_class import AssetClassManager
    from refgenie.managers.genome import GenomeManager


def provenance_columns(
    build_provenance: "BuildProvenance | None",
    peers: "Sequence[AssetName] | None" = None,
    group_label: str = "",
) -> dict[str, Any]:
    """
    ``AssetName`` provenance columns for a row about to be inserted.

    No provenance (the caller did not build) leaves every column NULL. That is
    the correct record for a name that no build produced, not a gap to fill.

    When the group already has a row for this exact ``build_digest``, the full
    provenance is still recorded on the new row: two names for one build is
    the same many-names-one-thing model the content layer already uses for
    ``Asset``, applied one level down. ``unique_group_build_digest`` is a
    non-unique index precisely so this insert does not collide.
    """
    if build_provenance is None:
        return {}
    if build_provenance.build_digest is not None and peers:
        already = next((p for p in peers if p.build_digest == build_provenance.build_digest), None)
        if already is not None:
            logger.info(
                f"Build {build_provenance.build_digest} is already recorded as "
                f"'{group_label}:{already.name}'; recording it again under this name too."
            )
    return build_provenance.as_columns()


def content_rel_path(genome_digest: GenomeDigest, asset_group_name: str, digest: str) -> Path:
    """The content-addressed data path, relative to the genome folder."""
    return Path("data") / genome_digest / asset_group_name / digest


def place_content(source: Path, content_abs: Path) -> bool:
    """
    Move freshly-produced content into its digest-addressed location.

    If the destination already exists, the content is a byte-for-byte
    duplicate that is already on disk, so the just-built/extracted ``source``
    is redundant and discarded. Colocation symlinks are relative and the
    source sits at the same depth as the destination, so they survive the
    move.

    The move is an ``os.rename``, which cannot cross filesystems. Source and
    destination are both under the (user-configurable) genome folder, so
    this holds unless something inside that folder is a mount point of its
    own -- in which case the rename raises ``OSError(EXDEV)`` rather than
    silently copying.

    Returns:
        bool: Whether content was moved, and can therefore be moved back.
    """
    if source.resolve() == content_abs.resolve():
        return False
    if content_abs.exists():
        if source.is_dir():
            shutil.rmtree(source)
        return False
    content_abs.parent.mkdir(parents=True, exist_ok=True)
    os.rename(source, content_abs)
    return True


def is_same_asset_group(asset: Asset, genome_digest: GenomeDigest, asset_group_name: str) -> bool:
    """Whether an asset belongs to the given genome's named asset group."""
    group = asset.asset_group
    if group is None:
        return False
    return group.name == asset_group_name and group.genome.digest == genome_digest


def exists_in_session(
    session: Session,
    genome_digest: GenomeDigest,
    asset_group_name: str,
    asset_name: str,
) -> bool:
    """
    Whether a name is taken in a group, asked of an already-open session.

    The same question as ``AssetManager.exists``, but for callers that already
    hold a session. ``_database_session`` constructs a new one on every access,
    so calling ``exists`` from inside an open block opens a second connection
    that cannot see the enclosing transaction and may contend with it.
    """
    return bool(
        session.exec(
            select(AssetName)
            .join(AssetGroup, AssetGroup.id == AssetName.asset_group_id)
            .join(Genome)
            .where(
                AssetName.name == asset_name,
                AssetGroup.name == asset_group_name,
                Genome.digest == genome_digest,
            )
        ).first()
    )


def ensure_asset_name(
    session: Session,
    asset: Asset,
    name: str,
    build_provenance: "BuildProvenance | None" = None,
) -> None:
    """Insert an AssetName row for (asset's group, name) unless it exists.

    Reached when a build completes an incomplete asset. The row may already
    exist -- ``add_incomplete`` writes one to keep the digest resolvable by
    name -- in which case it is the placeholder for this very build, and
    the provenance fills it in. A row that already records a build is left
    alone: that build is not this one.
    """
    peers = session.exec(
        select(AssetName).where(AssetName.asset_group_id == asset.asset_group_id)
    ).all()
    row = next((p for p in peers if p.name == name), None)
    columns = provenance_columns(
        build_provenance, peers=peers, group_label=str(asset.asset_group_id)
    )
    if row is None:
        session.add(
            AssetName(
                name=name,
                asset_group_id=asset.asset_group_id,
                asset_digest=asset.digest,
                **columns,
            )
        )
    elif columns and row.build_digest is None:
        for column, value in columns.items():
            setattr(row, column, value)


def named_asset_in_session(
    session: Session,
    genome_digest: GenomeDigest,
    asset_group_name: str,
    asset_name: str,
) -> Asset:
    """
    ``AssetManager.get``, asked of an already-open session.

    Raises:
        MissingAssetError: If no asset has that name.
    """
    return one_or_raise(
        session,
        asset_by_name_stmt(
            genome_digest,
            asset_group_name,
            asset_name,
            options=[
                selectinload(Asset.seek_keys),
                selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
                selectinload(Asset.asset_names),
            ],
        ),
        MissingAssetError(genome=genome_digest, asset_group=asset_group_name, asset=asset_name),
        unique=True,
    )


class AssetContentManager(ResourceManager):
    """
    The content write path: from produced content to catalog rows.

    ``add`` registers content (a build output, an extracted archive, a folder
    already under the genome folder), ``adopt_name`` gives content that is
    already here one more name, and ``add_incomplete`` fixes a digest before
    its content exists.

    Invariant: every name ``add`` or ``adopt_name`` commits has its alias tree
    rendered before the call returns; nothing renders before the commit.
    """

    def __init__(
        self,
        database_engine: Engine,
        genome_folder: Path,
        asset_class_manager: "AssetClassManager",
        genome_manager: "GenomeManager",
        groups: "AssetGroupManager",
        tree: "AliasTree",
        events: EventSink | None = None,
    ):
        """
        Initialize the AssetContentManager.

        Args:
            database_engine: The database engine.
            genome_folder: Path to genome data folder; content lands under its
                ``data/``.
            asset_class_manager: The AssetClassManager, for the class of added content.
            genome_manager: The GenomeManager, for the genome of a new group.
            groups: The AssetGroupManager, for group defaults.
            tree: The AliasTree that renders a newly named asset. ``add`` and
                ``adopt_name`` each render through it once, after their commit.
            events: Where committed asset changes are recorded for plugins.
        """
        super().__init__(database_engine)
        self._genome_folder = genome_folder
        self._asset_classes = asset_class_manager
        self._genomes = genome_manager
        self._groups = groups
        self._tree = tree
        self._events = events or NULL_EVENTS

    def _named_asset(
        self, genome_digest: GenomeDigest, asset_group_name: str, asset_name: str
    ) -> Asset:
        """The named asset, read in a session of its own."""
        with self._database_session as session:
            return named_asset_in_session(session, genome_digest, asset_group_name, asset_name)

    @update_scope
    def add(
        self,
        asset_class_name: str,
        path: Path,
        asset_group_name: str,
        *,
        genome_digest: GenomeDigest,
        asset_name: str | None = None,
        description: str | None = None,
        recipe: Recipe | None = None,
        custom_seek_keys: dict[str, "str | tuple[str, SeekKeyType]"] | None = None,
        colocate: list[dict[str, str]] | None = None,
        digest: str | None = None,
        set_default: bool | None = None,
        build_provenance: "BuildProvenance | None" = None,
    ) -> Asset | None:
        """
        Add an asset from a path.

        Args:
            asset_class_name: The name of the asset class.
            path: The path to the asset.
            asset_group_name: The name of the asset group.
            genome_digest: The genome digest.
            asset_name: The name of the asset.
            description: The description of the asset.
            recipe: The recipe used to build the asset.
            custom_seek_keys: Dict mapping seek key names to values. Values can be
                a plain string (defaults to SeekKeyType.string) or a
                (value, SeekKeyType) tuple for explicit type control.
            digest: The content digest, when it is already known and
                authoritative -- the pull path, where the server computed it
                under the building recipe's ``inherent`` set and the client
                cannot reproduce that. Omit it to compute the digest here.
            set_default: Whether to set this asset as the group default.
                True: always promote to default. False: never promote.
                None (default): promote only when creating a new group (the
                pull-path rule). The build path should pass True so deliberate builds
                become the default; the pull path should leave as None so
                pulled assets don't hijack the default.
            build_provenance: What the build that produced this content
                recorded about itself. Lands on the ``AssetName`` row, so a
                second build of identical content keeps its own. A caller that
                did not build passes nothing and the columns stay NULL.

        Returns:
            Asset: The added asset.
        """
        asset_class = self._asset_classes.get(asset_class_name)
        absolute_path = self._genome_folder / path
        if not absolute_path.exists():
            raise FileNotFoundError(
                f"Provided path does not exist: {absolute_path}. Please note that ",
                f"the provided path should be relative to the genome folder: {self._genome_folder}",
            )
        asset_name = asset_name or self._groups.get_default(
            asset_group_name, genome_digest=genome_digest
        )

        digest, content_rel, existing = self._resolve_target(
            absolute_path=absolute_path,
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            recipe=recipe,
            digest=digest,
        )

        # Every rejection decidable from metadata alone has now been made, so
        # the write is legal and the content may move. What remains can still
        # fail (a declared seek-key file may be missing, a concurrent writer may
        # win the insert), so the placement is undone if the catalog write does
        # not go through -- data/ has no sweeper, and content parked at a
        # digest path with no row would silently satisfy the next placement.
        content_abs = self._genome_folder / content_rel
        moved = place_content(absolute_path, content_abs)
        # Walked here rather than in an Asset before_insert listener: it is a
        # recursive glob over the whole asset, and running it inside the insert
        # flush held the write lock for its duration.
        content_size = directory_size(content_abs)
        try:
            if existing is not None:
                # The content is already in the catalog. Add a name for it (or
                # complete an incomplete row) rather than inserting a duplicate PK.
                asset = self._reconcile_existing_content(
                    existing=existing,
                    asset_name=asset_name,
                    content_rel=content_rel,
                    asset_class=asset_class,
                    custom_seek_keys=custom_seek_keys,
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    recipe=recipe,
                    description=description,
                    set_default=set_default,
                    content_size=content_size,
                    build_provenance=build_provenance,
                )
            else:
                asset = self._commit_new_content(
                    digest=digest,
                    content_rel=content_rel,
                    asset_class=asset_class,
                    custom_seek_keys=custom_seek_keys,
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    asset_name=asset_name,
                    recipe=recipe,
                    description=description,
                    colocate=colocate,
                    set_default=set_default,
                    content_size=content_size,
                    build_provenance=build_provenance,
                )
        except Exception:
            if moved:
                self._unplace_content(content_abs, absolute_path)
            raise
        # The name is committed; the alias tree is derived from it. Rendered
        # outside the try: a render failure must not unplace committed content.
        self._tree.render(genome_digest, asset_group_name, asset_name)
        return asset

    def _resolve_target(
        self,
        absolute_path: Path,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        recipe: Recipe | None = None,
        digest: str | None = None,
    ) -> tuple[str, Path, Asset | None]:
        """
        Identify the content and decide whether the catalog will accept it.

        Content addressing: the digest is the directory's identity, and the
        content belongs at data/<genome_digest>/<asset_group>/<content_digest>/.
        No canonical name appears anywhere in the content path.

        Only the files the recipe declares inherent contribute to the digest;
        incidental ones (logs, scratch) are still placed, just not identifying.
        getattr covers a duck-typed recipe; None there means every file counts.

        A caller that already holds an authoritative digest passes it instead.
        That is the pull path: the server computed the digest under the building
        recipe's ``inherent`` set, which does not travel with the archive, so
        recomputing here would disagree with the server for any recipe that
        declares one.

        This is also where every rejection that can be decided from the catalog
        alone is made, so that :meth:`add` can reject a write before
        moving anything on disk.

        Returns:
            tuple: (digest, content path relative to the genome folder, the
            catalog row already holding that digest or None).

        Raises:
            ValueError: If the digest belongs to another group, or the requested
                name in this group already maps to different content.
        """
        if digest is None:
            digest = get_dir_digest(absolute_path, inherent=getattr(recipe, "inherent", None))
        content_rel = content_rel_path(genome_digest, asset_group_name, digest)

        with self._database_session as session:
            existing = session.exec(asset_by_digest_stmt(digest)).unique().one_or_none()
            name_taken = existing is None and exists_in_session(
                session, genome_digest, asset_group_name, asset_name
            )
        if existing is not None:
            self._reject_cross_group(existing, genome_digest, asset_group_name)
        elif name_taken:
            # New content, but the requested name already maps to *different*
            # content in this group.
            raise ValueError(
                f"Asset name '{genome_digest}/{asset_group_name}:{asset_name}' already exists "
                f"and maps to different content."
            )
        return digest, content_rel, existing

    def _commit_new_content(
        self,
        digest: str,
        content_rel: Path,
        asset_class: Any,
        custom_seek_keys: dict[str, "str | tuple[str, SeekKeyType]"] | None,
        *,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        recipe: Recipe | None = None,
        description: str | None = None,
        colocate: list[dict[str, str]] | None = None,
        set_default: bool | None = None,
        content_size: int | None = None,
        build_provenance: "BuildProvenance | None" = None,
    ) -> Asset | None:
        """
        Insert a catalog row for content that is not yet known, and name it.

        The group row, the asset row, its seek keys and its name row are all
        written in one session, so a failure anywhere -- most commonly a seek
        key the asset class declares but the content does not satisfy -- leaves
        no partial catalog state behind.
        """
        is_new_group = not self._groups.exists(asset_group_name, genome_digest=genome_digest)
        # Resolve whether to promote this asset to group default:
        # - set_default=True: always promote (build path)
        # - set_default=False: never promote
        # - set_default=None: promote only for new groups (pull path)
        should_set_default = set_default if set_default is not None else is_new_group
        content_abs = self._genome_folder / content_rel

        # Resolved before the session opens: both of these read the database
        # through their own manager, and _database_session must not be
        # re-entered from inside an open block.
        if is_new_group:
            asset_group = AssetGroup(
                name=asset_group_name,
                genome=self._genomes.get(genome_digest),
                asset_class_id=asset_class.id,
                description=asset_class.description,
            )
        else:
            asset_group = self._groups.get(asset_group_name, genome_digest=genome_digest)

        with self._database_session as session:
            if is_new_group:
                session.add(asset_group)
            asset = Asset(
                name=asset_name,
                asset_group=asset_group,
                path=content_rel.as_posix(),
                digest=digest,
                description=description,
                recipe_id=recipe.id if recipe else None,
                colocate=colocate,
                size=content_size,
            )
            session.add(asset)
            bind_asset_class_seek_keys(
                asset=asset,
                asset_class=asset_class,
                directory_path=content_abs,
                custom_seek_keys=custom_seek_keys,
                session=session,
            )
            session.flush()
            # The name axis: a row mapping (asset_group, name) -> content
            # digest, carrying the provenance of the build that produced it.
            session.add(
                AssetName(
                    name=asset_name,
                    asset_group_id=asset.asset_group_id,
                    asset_digest=asset.digest,
                    **provenance_columns(build_provenance),
                )
            )
            try:
                session.commit()
            except IntegrityError:
                # Concurrency backstop: with no locking anywhere in refgenie, two
                # different targets producing byte-identical content concurrently
                # can both reach here. The digest PK / name uniqueness rejects the
                # second insert; re-read and return the winner.
                session.rollback()
                logger.info(
                    f"Concurrent insert for '{genome_digest}/{asset_group_name}:{asset_name}'; "
                    f"returning the existing row."
                )
                return named_asset_in_session(session, genome_digest, asset_group_name, asset_name)
            session.refresh(asset)
            self._record_added(genome_digest, asset_group_name, asset_name, asset.digest)
            # Eagerly load the asset with its relationships before returning
            asset = (
                session.exec(
                    select(Asset)
                    .where(Asset.digest == asset.digest)
                    .options(
                        selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
                        selectinload(Asset.asset_names),
                    )
                )
                .unique()
                .one()
            )
        logger.info(f"Added: '{genome_digest}/{asset_group_name}:{asset_name}'")
        if should_set_default:
            self._groups.set_default(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
            )
        return asset

    def _unplace_content(self, content_abs: Path, source: Path) -> None:
        """
        Undo a ``place_content`` move after the catalog write failed.

        The directories the placement created on the way in are pruned too, so a
        rejected write leaves no ``data/<genome>/<group>/`` for a group that was
        never created. Pruning stops at the first non-empty directory.

        Best effort: a failure to restore is logged rather than raised, because
        it must not mask the exception that triggered the rollback.
        """
        try:
            os.rename(content_abs, source)
        except OSError as e:
            logger.error(
                f"Could not restore {content_abs} to {source} after a failed add: {e}. "
                f"The content is orphaned at its digest path with no catalog row."
            )
            return
        parent = content_abs.parent
        data_folder = self._genome_folder / "data"
        while parent != data_folder and parent.is_relative_to(data_folder):
            try:
                parent.rmdir()
            except OSError:
                break
            parent = parent.parent

    def _reconcile_existing_content(
        self,
        existing: Asset,
        asset_name: str,
        content_rel: Path,
        asset_class: Any,
        custom_seek_keys: dict[str, "str | tuple[str, SeekKeyType]"] | None,
        *,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        recipe: Any = None,
        description: str | None = None,
        set_default: bool | None = None,
        content_size: int | None = None,
        build_provenance: "BuildProvenance | None" = None,
    ) -> Asset:
        """
        Attach a name to content that is already in the catalog.

        Two cases (the third, a cross-group collision, is rejected upstream by
        :meth:`_resolve_target` before the content moves; the guard is repeated
        here so the method is safe to call on its own):
          * the existing row is *incomplete* (digest set in stone, no path) -- it
            is completed by setting its path and binding seek keys;
          * the existing row is complete -- the requested name is added as a peer
            ``AssetName`` (or is a no-op if the name already points here).
        """
        self._reject_cross_group(existing, genome_digest, asset_group_name)

        content_abs = self._genome_folder / content_rel

        if existing.path is None:
            # Complete an incomplete asset.
            with self._database_session as session:
                asset = (
                    session.exec(select(Asset).where(Asset.digest == existing.digest))
                    .unique()
                    .one()
                )
                asset.path = content_rel.as_posix()
                # Completing an incomplete asset is an UPDATE, so no insert
                # listener ever ran for it and its size was left null.
                asset.size = content_size
                if recipe is not None:
                    asset.recipe_id = recipe.id
                if description is not None:
                    asset.description = description
                bind_asset_class_seek_keys(
                    asset=asset,
                    asset_class=asset_class,
                    directory_path=content_abs,
                    custom_seek_keys=custom_seek_keys,
                    session=session,
                )
                ensure_asset_name(session, asset, asset_name, build_provenance)
                session.commit()
            logger.info(
                f"Completed incomplete asset '{genome_digest}/{asset_group_name}:{asset_name}'"
            )
            self._record_added(genome_digest, asset_group_name, asset_name, existing.digest)
            if set_default is True:
                self._groups.set_default(
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    asset_name=asset_name,
                )
            return self._named_asset(genome_digest, asset_group_name, asset_name)

        # Complete row already holds this content: add the requested name if new.
        with self._database_session as session:
            asset = (
                session.exec(
                    select(Asset)
                    .where(Asset.digest == existing.digest)
                    .options(selectinload(Asset.asset_names))
                )
                .unique()
                .one()
            )
            if any(n.name == asset_name for n in asset.asset_names):
                logger.info(
                    f"Asset '{genome_digest}/{asset_group_name}:{asset_name}' already exists with "
                    f"identical content (digest {asset.digest}); nothing to do."
                )
            else:
                # Reuse the open session: _database_session builds a new one per
                # access, and a second connection could neither see this
                # transaction nor safely contend with it.
                if exists_in_session(
                    session=session,
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    asset_name=asset_name,
                ):
                    raise ValueError(
                        f"Asset name '{genome_digest}/{asset_group_name}:{asset_name}' already "
                        f"exists and maps to different content."
                    )
                # A second build of identical content under a new name: the
                # content row already exists, so this build's provenance must
                # go on the new name row or it is lost.
                session.add(
                    AssetName(
                        name=asset_name,
                        asset_group_id=asset.asset_group_id,
                        asset_digest=asset.digest,
                        **provenance_columns(
                            build_provenance,
                            peers=asset.asset_names,
                            group_label=f"{genome_digest}/{asset_group_name}",
                        ),
                    )
                )
                session.commit()
                logger.info(
                    f"Added name '{genome_digest}/{asset_group_name}:{asset_name}' -> "
                    f"existing content {asset.digest}"
                )
                self._record_added(genome_digest, asset_group_name, asset_name, asset.digest)
        if set_default is True:
            self._groups.set_default(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
            )
        return self._named_asset(genome_digest, asset_group_name, asset_name)

    @update_scope
    def adopt_name(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        asset_digest: str,
        build_provenance: "BuildProvenance | None" = None,
    ) -> None:
        """
        Attach an additional name to content already present locally.

        Used by pull to record every name a server reports for one content
        digest, along with the build the server recorded for that name -- so a
        pulled build is addressable by its ``build_digest`` on the client too,
        not only on the machine that built it. A name already taken by
        *different* content in the group is skipped with a warning rather than
        raising, so a name clash never fails a pull.
        """
        with self._database_session as session:
            existing = session.exec(
                select(AssetName)
                .join(AssetGroup, AssetGroup.id == AssetName.asset_group_id)
                .join(Genome)
                .where(
                    AssetName.name == asset_name,
                    AssetGroup.name == asset_group_name,
                    Genome.digest == genome_digest,
                )
            ).one_or_none()
            if existing is not None:
                if existing.asset_digest != asset_digest:
                    logger.warning(
                        f"Name '{genome_digest}/{asset_group_name}:{asset_name}' already maps to "
                        f"different content ({existing.asset_digest}); not adopting."
                    )
                return
            asset = (
                session.exec(select(Asset).where(Asset.digest == asset_digest))
                .unique()
                .one_or_none()
            )
            if asset is None:
                logger.warning(f"Cannot adopt name for missing content {asset_digest}")
                return
            peers = session.exec(
                select(AssetName).where(AssetName.asset_group_id == asset.asset_group_id)
            ).all()
            session.add(
                AssetName(
                    name=asset_name,
                    asset_group_id=asset.asset_group_id,
                    asset_digest=asset_digest,
                    **provenance_columns(
                        build_provenance,
                        peers=peers,
                        group_label=f"{genome_digest}/{asset_group_name}",
                    ),
                )
            )
            session.commit()
        self._record_added(genome_digest, asset_group_name, asset_name, asset_digest)
        self._tree.render(genome_digest, asset_group_name, asset_name)

    def _record_added(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        digest: str,
    ) -> None:
        """Record that ``asset_name`` now resolves to committed content."""
        self._events.record(
            Change(
                action="asset_added",
                genome=genome_digest,
                asset_group=asset_group_name,
                asset=asset_name,
                digest=digest,
            )
        )

    def _reject_cross_group(
        self,
        existing: Asset,
        genome_digest: GenomeDigest,
        asset_group_name: str,
    ) -> None:
        """Raise if content already in the catalog belongs to a different group."""
        if not is_same_asset_group(existing, genome_digest, asset_group_name):
            raise ValueError(
                f"Asset digest {existing.digest} is already used by '{existing.registry_path}', "
                f"which is not in '{genome_digest}/{asset_group_name}'. Cross-group content "
                f"collisions are not representable while asset.digest is a global primary key."
            )

    def add_incomplete(
        self,
        digest: str,
        asset_class_name: str,
        asset_group_name: str,
        *,
        genome_digest: GenomeDigest,
        asset_name: str | None = None,
        description: str | None = None,
    ) -> Asset | None:
        """
        Add an incomplete asset (no path attribute defined).

        Incomplete assets are used to set the digest in stone and prevent
        breaking children-parent provenance chain.

        Args:
            digest: The digest of the asset.
            asset_class_name: The name of the asset class.
            asset_group_name: The name of the asset group.
            genome_digest: The genome digest.
            asset_name: The name of the asset.
            description: The description of the asset.

        Returns:
            Asset: The added incomplete asset.
        """
        asset_class = self._asset_classes.get(asset_class_name)
        asset_name = asset_name or self._groups.get_default(
            asset_group_name, genome_digest=genome_digest
        )
        with self._database_session as session:
            if exists_in_session(session, genome_digest, asset_group_name, asset_name):
                raise ValueError(
                    f"Asset '{genome_digest}/{asset_group_name}:{asset_name}' already exists"
                )
            try:
                asset_group = asset_group_in_session(
                    session, asset_group_name, genome_digest=genome_digest
                )
            except MissingAssetGroupError:
                asset_group = AssetGroup(
                    name=asset_group_name,
                    genome=one_or_raise(
                        session,
                        select(Genome).where(Genome.digest == genome_digest),
                        MissingGenomeError(genome=genome_digest),
                    ),
                    asset_class=asset_class,
                )
                session.add(asset_group)
            asset = Asset(
                name=asset_name,
                asset_group=asset_group,
                digest=digest,
                description=description,
            )
            session.add(asset)
            session.flush()
            # Even an incomplete asset needs a name row, or it is unresolvable
            # by name (get/exists/seek all route through assetname).
            session.add(
                AssetName(
                    name=asset_name,
                    asset_group_id=asset.asset_group_id,
                    asset_digest=asset.digest,
                )
            )
            session.commit()
            # Inside the session: the name lazy-loads the group's genome.
            logger.info(f"Added incomplete asset: {asset}")
        self._groups.set_default(asset_group_name, asset_name, genome_digest=genome_digest)
        return asset
