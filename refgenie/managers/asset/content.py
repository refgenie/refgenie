"""
AssetContentManager - the content write path, plus its stateless helpers.

Reached as ``rgc.asset.content``. It turns produced content (a build output, a
pulled archive, a folder the user registers) into catalog rows: content
addressing, placement under ``data/``, and the name rows. Its reads go straight
to ``refgenie.db.queries``, so it holds no reference to ``AssetManager``.

The module functions run inside a caller's open session or use no manager
state at all.
"""

import os
import shutil
from pathlib import Path
from typing import Any, TYPE_CHECKING

from refgenie.db.engine import Engine
from refgenie.db.queries import asset as queries
from refgenie.db.session import Session, commit_or_conflict
from refgenie.db.tables import Asset, AssetGroup, AssetName, Recipe, SeekKeyType
from refgenie.exceptions import CatalogConflictError
from refgenie.logger import logger
from refgenie.managers.asset.seek_key import bind_asset_class_seek_keys
from refgenie.managers.base import ResourceManager
from refgenie.models import GenomeDigest
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import Change
from refgenie.utils.digest import BuildProvenance, directory_size, get_dir_digest

if TYPE_CHECKING:
    from refgenie.managers.asset.alias_tree import AliasTree
    from refgenie.managers.asset.group import AssetGroupManager
    from refgenie.managers.asset_class import AssetClassManager
    from refgenie.managers.genome import GenomeManager


def content_rel_path(genome_digest: GenomeDigest, asset_group_name: str, digest: str) -> Path:
    """The content-addressed data path, relative to the genome folder."""
    return Path("data") / genome_digest / asset_group_name / digest


def place_content(source: Path, content_abs: Path) -> bool:
    """
    Move freshly-produced content into its digest-addressed location.

    If the destination already exists, it holds identical content, so
    ``source`` is discarded. The move is an ``os.rename``, so it raises
    ``OSError(EXDEV)`` if the genome folder spans filesystems.

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


def _add_name(
    session: Session, asset: Asset, name: str, build_provenance: "BuildProvenance | None"
) -> None:
    """Insert the ``AssetName`` row giving ``asset`` the name ``name``."""
    session.add(
        AssetName(
            name=name,
            asset_group_id=asset.asset_group_id,
            asset_digest=asset.digest,
            **(build_provenance.as_columns() if build_provenance else {}),
        )
    )


def _name_taken(genome_digest: GenomeDigest, asset_group_name: str, asset_name: str) -> ValueError:
    return ValueError(
        f"Asset name '{genome_digest}/{asset_group_name}:{asset_name}' already exists "
        f"and maps to different content."
    )


class AssetContentManager(ResourceManager):
    """
    The content write path: from produced content to catalog rows.

    ``add`` registers content (a build output, an extracted archive, a folder
    already under the genome folder), and ``adopt_name`` gives content that is
    already here one more name.

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
            return queries.get_by_name(session, genome_digest, asset_group_name, asset_name)

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
            digest: The content digest, when already known (pull passes the
                server's, which the client cannot recompute). Omit to compute it.
            set_default: Whether to make this asset the group default. True:
                always (build). False: never. None: only for a new group (pull).
            build_provenance: The producing build's record, stored on the
                ``AssetName`` row. Omit when nothing was built.

        Returns:
            Asset: The added asset.
        """
        asset_class = self._asset_classes.get(asset_class_name)
        absolute_path = self._genome_folder / path
        if not absolute_path.exists():
            raise FileNotFoundError(
                f"Provided path does not exist: {absolute_path}. "
                f"The path must be relative to the genome folder: {self._genome_folder}"
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

        # The catalog write can still fail, so the placement is undone on error:
        # orphaned content at a digest path would satisfy the next placement.
        content_abs = self._genome_folder / content_rel
        moved = place_content(absolute_path, content_abs)
        # Measured outside the insert flush so the walk does not hold the write lock.
        content_size = directory_size(content_abs)
        try:
            if existing is not None:
                asset = self._reconcile_existing_content(
                    existing=existing,
                    asset_name=asset_name,
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    set_default=set_default,
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
        # Outside the try: a render failure must not unplace committed content.
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

        The digest covers only the recipe's ``inherent`` files (all files if it
        declares none), and the content belongs at
        ``data/<genome_digest>/<asset_group>/<digest>/``. Every rejection the
        catalog alone can decide is made here, before anything moves on disk.

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
            existing = queries.find_by_digest(session, digest)
            name_taken = existing is None and queries.exists_name(
                session, genome_digest, asset_group_name, asset_name
            )
        if existing is not None:
            self._reject_cross_group(existing, genome_digest, asset_group_name)
        elif name_taken:
            raise _name_taken(genome_digest, asset_group_name, asset_name)
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

        All rows are written in one session, so a failure (often a missing
        seek-key file) leaves no partial catalog state.
        """
        is_new_group = not self._groups.exists(asset_group_name, genome_digest=genome_digest)
        should_set_default = set_default if set_default is not None else is_new_group
        content_abs = self._genome_folder / content_rel

        # Read before the session opens: _database_session must not be re-entered.
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
            _add_name(session, asset, asset_name, build_provenance)
            try:
                commit_or_conflict(session)
            except CatalogConflictError:
                # Two concurrent writers of identical content: the loser re-reads the winner.
                logger.info(
                    f"Concurrent insert for '{genome_digest}/{asset_group_name}:{asset_name}'; "
                    f"returning the existing row."
                )
                return queries.get_by_name(session, genome_digest, asset_group_name, asset_name)
            session.refresh(asset)
            self._record_added(genome_digest, asset_group_name, asset_name, asset.digest)
            asset = queries.get_by_digest(session, asset.digest)
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

        Also prunes the now-empty directories the placement created. A failure
        to restore is logged, not raised, so it does not mask the original error.
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
        *,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        set_default: bool | None = None,
        build_provenance: "BuildProvenance | None" = None,
    ) -> Asset:
        """Give content already in the catalog this name, unless it has it already."""
        with self._database_session as session:
            asset = queries.get_by_digest(session, existing.digest)
            if any(n.name == asset_name for n in asset.asset_names):
                logger.info(
                    f"Asset '{genome_digest}/{asset_group_name}:{asset_name}' already exists with "
                    f"identical content (digest {asset.digest}); nothing to do."
                )
            else:
                if queries.exists_name(session, genome_digest, asset_group_name, asset_name):
                    raise _name_taken(genome_digest, asset_group_name, asset_name)
                _add_name(session, asset, asset_name, build_provenance)
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

        Pull uses this to record every name a server reports for a digest, with
        its build provenance. A name taken by different content is skipped with
        a warning, so a clash never fails a pull.
        """
        with self._database_session as session:
            existing = queries.find_name(session, genome_digest, asset_group_name, asset_name)
            if existing is not None:
                if existing.asset_digest != asset_digest:
                    logger.warning(
                        f"Name '{genome_digest}/{asset_group_name}:{asset_name}' already maps to "
                        f"different content ({existing.asset_digest}); not adopting."
                    )
                return
            asset = queries.find_by_digest(session, asset_digest)
            if asset is None:
                logger.warning(f"Cannot adopt name for missing content {asset_digest}")
                return
            _add_name(session, asset, asset_name, build_provenance)
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
