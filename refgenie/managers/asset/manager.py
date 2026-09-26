"""
AssetManager - asset records, listing, removal, rename and local seek.

AssetContentManager, AssetGroupManager, SeekKeyManager, AliasTree and
AssetLinkManager are composed in as the public attributes ``content``,
``group``, ``seek_key``, ``tree`` and ``links``.
Building and pulling are not here: ``Refgenie`` owns the builder and the
puller, and both depend on this manager, never the other way round.
"""

import shutil
from pathlib import Path
from collections.abc import Iterable
from typing import TYPE_CHECKING
from collections import defaultdict

from rich.table import Table
from sqlalchemy.engine import Engine
from sqlalchemy.orm import selectinload
from sqlmodel import select

from refgenie.db.tables import Asset, AssetGroup, AssetName, Genome, is_path_type
from refgenie.exceptions import MissingAssetError, MissingAssetGroupError
from refgenie.logger import logger
from refgenie.models import AssetRegistryPathComponents, GenomeAlias, GenomeDigest
from refgenie.managers.base import ResourceManager
from refgenie.managers.asset.alias_tree import AliasTree
from refgenie.managers.asset.content import (
    AssetContentManager,
    exists_in_session,
    named_asset_in_session,
)
from refgenie.managers.asset.group import AssetGroupManager
from refgenie.managers.asset.links import AssetLinkManager
from refgenie.managers.asset.queries import asset_by_digest_stmt, asset_by_name_stmt, assets_stmt
from refgenie.managers.asset.seek_key import SeekKeyManager
from refgenie.managers.asset.tables import asset_table
from refgenie.managers.queries import one_or_raise
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import Change

if TYPE_CHECKING:
    from refgenie.managers.alias import AliasManager
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.asset_class import AssetClassManager


class AssetManager(ResourceManager):
    """
    Manager for asset records: lookups, listing, removal, rename and local seek.

    Five composed managers are public attributes, used directly rather than
    through wrapper methods: ``content`` (the content write path: ``add``,
    ``adopt_name``, ``add_incomplete``), ``group`` (asset groups and their
    default asset), ``seek_key`` (seek-key lookups), ``tree`` (the per-alias
    ``alias/`` and ``builds/`` trees on disk) and ``links`` (parent/child links).

    Remote listing and remote seek are ``rgc.servers``; building and pulling
    are ``rgc.build.run`` / ``rgc.transfer.pull``.
    """

    def __init__(
        self,
        database_engine: Engine,
        genome_folder: Path,
        alias_folder: Path,
        alias_manager: "AliasManager",
        genome_manager: "GenomeManager",
        asset_class_manager: "AssetClassManager",
        events: EventSink | None = None,
    ):
        """
        Initialize the AssetManager.

        Args:
            database_engine: The database engine.
            genome_folder: Path to genome data folder.
            alias_folder: Path to alias folder.
            alias_manager: The AliasManager for resolving genome names.
            genome_manager: The GenomeManager, for ``group`` and ``content``.
            asset_class_manager: The AssetClassManager, for ``content``.
            events: Where committed asset changes are recorded for plugins.
                Shared with ``group`` and ``content``.
        """
        super().__init__(database_engine)
        self._events = events or NULL_EVENTS
        self._genome_folder = genome_folder
        self._alias_folder = alias_folder
        self._alias_manager = alias_manager
        self.group = AssetGroupManager(database_engine, genome_manager, events=self._events)
        self.seek_key = SeekKeyManager(database_engine, self.group)
        self.tree = AliasTree(
            database_engine, genome_folder, alias_folder, alias_manager, self.group
        )
        self.links = AssetLinkManager(database_engine)
        self.content = AssetContentManager(
            database_engine,
            genome_folder,
            asset_class_manager,
            genome_manager,
            self.group,
            self.tree,
            events=self._events,
        )

    @property
    def genome_folder(self) -> Path:
        """Get the genome folder path."""
        return self._genome_folder

    @property
    def data_folder(self) -> Path:
        """Get the data directory."""
        return self._genome_folder / "data"

    @property
    def alias_folder(self) -> Path:
        """Get the alias directory."""
        return self._alias_folder

    # === CRUD Operations ===

    def get_by_digest(self, digest: str) -> Asset | None:
        """
        Get an asset by its content digest, or None if no asset holds that digest.

        Returns None rather than raising when there is no match, so it also
        answers whether a piece of content is already in the catalog.

        Args:
            digest: The digest of the asset.

        Returns:
            Asset | None: The asset, or None if not found.
        """
        with self._database_session as session:
            return session.exec(asset_by_digest_stmt(digest)).unique().one_or_none()

    def get(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> Asset:
        """
        Get an asset by its name. For content by digest, use :meth:`get_by_digest`.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The genome digest.

        Returns:
            Asset: The asset.

        Raises:
            MissingAssetError: If no asset has that name.
        """
        with self._database_session as session:
            return named_asset_in_session(session, genome_digest, asset_group_name, asset_name)

    def get_asset_dir(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> Path:
        """
        Absolute path to the asset's data directory.

        This is the single definition of where an asset's files live. Callers
        that need the directory must use this rather than deriving it from a
        seek key, whose value may be nested and whose type may not be a path.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The genome digest.

        Returns:
            Path: The absolute asset directory.
        """
        asset = self.get(
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            genome_digest=genome_digest,
        )
        if asset.path is None:
            raise ValueError(f"Incomplete asset, path is not set: {asset.registry_path}")
        return self._genome_folder / asset.path

    @update_scope
    def remove(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_digest: GenomeDigest,
        keep_asset_group: bool = False,
    ) -> str:
        """
        Remove an asset's metadata and files, addressed by name.

        This is the name axis onto :meth:`remove_by_digest`: it resolves the name
        to a content digest and delegates. Both entry points therefore perform
        exactly one removal algorithm.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The genome digest.
            keep_asset_group: Whether to keep the asset group if no assets left.

        Returns:
            The registry path of the removed asset.

        Raises:
            ValueError: If ``asset_name`` is not the content's canonical name.
            MissingAssetError: If the asset does not exist.
        """

        with self._database_session as session:
            # Resolve by name through the assetname table.
            asset = one_or_raise(
                session,
                asset_by_name_stmt(
                    genome_digest,
                    asset_group_name,
                    asset_name,
                    options=[selectinload(Asset.asset_names)],
                ),
                MissingAssetError(
                    genome=genome_digest, asset_group=asset_group_name, asset=asset_name
                ),
                unique=True,
            )

            # Removal destroys the content and *every* name that points at it.
            # Require the canonical (publication) name, and report what would be
            # lost otherwise, so a caller cannot silently wipe peer names via a
            # secondary one.
            if asset_name != asset.name:
                asset_names = [n.name for n in asset.asset_names]
                raise ValueError(
                    f"Can't remove via non-canonical name "
                    f"'{genome_digest}/{asset_group_name}:{asset_name}'. This content is known "
                    f"as {asset_names}; removing it would destroy all of them. Use the canonical "
                    f"name '{asset.name}'."
                )
            digest = asset.digest

        return self.remove_by_digest(digest=digest, keep_asset_group=keep_asset_group)

    @update_scope
    def remove_by_digest(self, digest: str, keep_asset_group: bool = False) -> str:
        """
        Remove an asset by its digest, including symlink cleanup and group removal.

        This is the single implementation of asset removal; :meth:`remove`
        resolves a name to a digest and calls it.

        Args:
            digest: The asset digest.
            keep_asset_group: Whether to keep the asset group if no assets left.

        Returns:
            The registry path of the removed asset.

        Raises:
            ValueError: If the asset has children.
            MissingAssetError: If the asset does not exist.
        """
        with self._database_session as session:
            asset = one_or_raise(
                session,
                select(Asset)
                .where(Asset.digest == digest)
                .options(
                    selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
                    selectinload(Asset.children),
                    selectinload(Asset.asset_names),
                ),
                MissingAssetError(digest=digest),
                unique=True,
            )

            if asset.children:
                raise ValueError(
                    f"Can't remove. Asset '{asset.registry_path}' has children. "
                    f"Remove them first: {', '.join(c.registry_path for c in asset.children)}"
                )

            # Extract values before deletion
            genome_digest = asset.asset_group.genome.digest
            asset_group_name = asset.asset_group.name
            asset_names = [n.name for n in asset.asset_names]
            asset_digest = asset.digest
            registry_path = asset.registry_path

            # Get asset group for checking emptiness
            asset_group = one_or_raise(
                session,
                select(AssetGroup)
                .options(selectinload(AssetGroup.assets))
                .join(Genome)
                .where(AssetGroup.name == asset_group_name, Genome.digest == genome_digest),
                MissingAssetGroupError(genome=genome_digest, asset_group=asset_group_name),
                unique=True,
            )

            canonical_name = asset.name
            session.delete(asset)
            session.commit()

            remaining = [a for a in asset_group.assets if a.digest != asset_digest]

        logger.info(f"Removed asset '{registry_path}'")
        self._events.record(
            Change(
                action="asset_removed",
                genome=genome_digest,
                asset_group=asset_group_name,
                asset=canonical_name,
                digest=asset_digest,
            )
        )
        self.tree.purge_names(genome_digest, asset_group_name, asset_names)

        if remaining or keep_asset_group:
            return registry_path

        # Empty group cleanup. The group's directories are resolved *before* the
        # group row goes away: removing the last group of a genome removes the
        # genome, and with it the aliases these paths are keyed by.
        logger.info(f"No assets left for '{genome_digest}/{asset_group_name}'. Removing group.")
        group_trees = self.tree.owned_trees(genome_digest, asset_group_name)
        self.group.remove_rows(asset_group_name, genome_digest)
        for alias_name, path in group_trees:
            logger.info(f"Removing group files for alias '{alias_name}': {path}")
            shutil.rmtree(path, ignore_errors=True)

        return registry_path

    def exists(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> bool:
        """
        Check if an asset name exists. For content by digest, use :meth:`get_by_digest`.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The genome digest.

        Returns:
            bool: Whether the asset exists. False for an unknown genome.
        """
        with self._database_session as session:
            return exists_in_session(session, genome_digest, asset_group_name, asset_name)

    def list_assets(
        self,
        genome_digests: list[GenomeDigest] | None = None,
        asset_group_name: str | None = None,
    ) -> Iterable[Asset]:
        """
        List all assets.

        Args:
            genome_digests: Genome digests to filter by.
            asset_group_name: The name of asset group to filter by.

        Returns:
            Iterable[Asset]: A list of all assets.
        """
        statement = assets_stmt(genome_digests)
        if asset_group_name is not None:
            statement = statement.where(AssetGroup.name == asset_group_name)

        with self._database_session as session:
            result = session.exec(statement)
            return result.unique().all()

    def list_seek_keys_values(
        self,
        genome_digests: list[GenomeDigest] | None = None,
        asset_group_name: str | None = None,
    ) -> dict[str, dict[str, dict[str, dict[str, str]]]]:
        """
        Bulk-load all seek key values, returning a nested mapping suitable
        for populator-style consumers.

        The shape is:
            {genome_name: {asset_group_name: {asset_name: {seek_key_name: str}}}}

        For path-type seek keys the value is the absolute path string;
        non-path seek keys return their stored value as a string. All values
        are guaranteed to be ``str`` (no ``Path`` instances) so the result is
        directly JSON-serializable and Jinja-templatable.

        Args:
            genome_digests: Optional genome digests to filter by.
            asset_group_name: Optional asset group name to filter by.

        Returns:
            A nested dict keyed by primary genome alias, asset group, asset,
            and seek key name with str values.
        """
        assets = self.list_assets(genome_digests=genome_digests, asset_group_name=asset_group_name)

        # Map genome_digest -> primary alias name (first alias) once
        digest_to_alias: dict[str, str] = {}

        result: dict[str, dict[str, dict[str, dict[str, str]]]] = {}
        for asset in assets:
            asset_group = asset.asset_group
            genome = asset_group.genome
            genome_digest = genome.digest
            if genome_digest not in digest_to_alias:
                aliases = self._alias_manager.get_for_genome(genome_digest)
                digest_to_alias[genome_digest] = aliases[0] if aliases else genome_digest
            genome_key = digest_to_alias[genome_digest]

            per_asset: dict[str, str] = {}
            for sk in asset.seek_keys:
                if is_path_type(sk.type):
                    if asset.path is None:
                        continue
                    per_asset[sk.name] = str(
                        (self._genome_folder / asset.path / sk.value).absolute()
                    )
                else:
                    per_asset[sk.name] = str(sk.value)

            result.setdefault(genome_key, {}).setdefault(asset_group.name, {})[asset.name] = (
                per_asset
            )

        return result

    def list_all(
        self,
        genome_digests: list[GenomeDigest] | None = None,
        asset_group_name: str | None = None,
        include_seek_keys: bool = False,
    ) -> tuple[dict[str, list[str]], dict[str, str]]:
        """
        List assets with formatted output.

        Args:
            genome_digests: Genome digests to filter by.
            asset_group_name: The name of the asset group.
            include_seek_keys: Whether to include seek keys.

        Returns:
            tuple[dict[str, list[str]], dict[str, str]]: Asset data and aliases data.
        """
        assets = self.list_assets(genome_digests=genome_digests, asset_group_name=asset_group_name)
        return_dict = defaultdict(list)
        for asset in assets:
            asset_strings = (
                [
                    f"{asset.asset_group.name}.{seek_key.name}:{asset.name}"
                    for seek_key in asset.seek_keys
                ]
                if include_seek_keys
                else [f"{asset.asset_group.name}:{asset.name}"]
            )
            return_dict[asset.asset_group.genome.digest].extend(asset_strings)

        # Get aliases for the genomes that have assets
        genome_digests_with_assets = list(return_dict.keys())
        aliases_dict = {}
        if genome_digests_with_assets:
            local_aliases = self._alias_manager.list_all()
            # Group aliases by genome digest
            aliases_by_genome = defaultdict(list)
            for alias in local_aliases:
                if alias.genome_digest in genome_digests_with_assets:
                    aliases_by_genome[alias.genome_digest].append(alias.name)

            # Convert to comma-separated strings
            for genome_digest, alias_names in aliases_by_genome.items():
                aliases_dict[genome_digest] = ", ".join(alias_names)

        return dict(return_dict), aliases_dict

    def table(
        self,
        genome_digests: list[GenomeDigest] | None = None,
        include_seek_keys: bool = False,
    ) -> list[Table]:
        """
        Create a table of all assets.

        Args:
            genome_digests: Genome digests to filter by.
            include_seek_keys: Whether to include seek keys.

        Returns:
            list[Table]: A list of tables.
        """
        asset_data, aliases_data = self.list_all(
            genome_digests=genome_digests, include_seek_keys=include_seek_keys
        )
        return [
            asset_table(
                asset_data=asset_data,
                aliases_data=aliases_data,
                include_seek_keys=include_seek_keys,
                source="local",
            )
        ]

    @update_scope
    def rename(
        self,
        asset_group_name: str,
        asset_name: str,
        new_asset_name: str,
        *,
        genome_digest: GenomeDigest,
    ):
        """
        Rename an asset.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            new_asset_name: The new name for the asset.
            genome_digest: The genome digest.

        Returns:
            The renamed asset object.
        """
        asset_group = self.group.get(asset_group_name, genome_digest=genome_digest)

        with self._database_session as session:
            row = one_or_raise(
                session,
                select(AssetName).where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.name == asset_name,
                ),
                MissingAssetError(
                    genome=genome_digest, asset_group=asset_group_name, asset=asset_name
                ),
            )
            if session.exec(
                select(AssetName).where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.name == new_asset_name,
                )
            ).first():
                raise ValueError(
                    f"Asset name '{genome_digest}/{asset_group_name}:{new_asset_name}' "
                    f"already exists."
                )
            asset = one_or_raise(
                session,
                select(Asset).where(Asset.digest == row.asset_digest),
                MissingAssetError(digest=row.asset_digest),
                unique=True,
            )
            # Renaming the name axis. Content does not move (it is
            # digest-addressed). If this is the publication name, carry it too so
            # the staged-tarball/S3 key follows.
            was_canonical = asset.name == asset_name
            row.name = new_asset_name
            if was_canonical:
                asset.name = new_asset_name
            session.add(row)
            session.add(asset)
            session.commit()
            asset_digest = asset.digest
        self._events.record(
            Change(
                action="asset_renamed",
                genome=genome_digest,
                asset_group=asset_group_name,
                asset=new_asset_name,
                previous=asset_name,
                digest=asset_digest,
            )
        )

        # Render the alias tree under the new name and clear the stale trees the
        # old name owned. Both of them: the build directory holds the completion
        # flag, and leaving it behind means a later build of the OLD name is
        # skipped as already done, while a rename back to it fails with "already
        # exists". owned_trees is the same enumeration purge_names uses,
        # so the two removal paths cannot disagree about what a name owns.
        self.tree.render(genome_digest, asset_group_name, new_asset_name)
        for alias_name, path in self.tree.owned_trees(genome_digest, asset_group_name, asset_name):
            logger.info(f"Removing stale files for alias '{alias_name}': {path}")
            shutil.rmtree(path, ignore_errors=True)
        logger.info(
            f"Renamed asset '{genome_digest}/{asset_group_name}:{asset_name}' to '{new_asset_name}'"
        )
        return self.get(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=new_asset_name,
        )

    # =========================================================================
    # Seek: registry path -> local alias-tree/content path, or remote file URL
    # =========================================================================

    def seek(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None = None,
        seek_key_name: str | None = None,
        force_exists: bool = False,
        abs_path: bool = False,
    ) -> str:
        """
        Get the path to an asset file, or the value of a non-path seek key.

        By default the returned path is under the human-readable alias tree
        (``alias/<alias>/<group>/<name>/...``), named for the genome's first
        local alias; a genome with no alias gets the content path. Pass
        ``abs_path=True`` for the digest-addressed content path under ``data/``
        instead. To name the alias the path is keyed by, use
        :meth:`seek_components` with a registry path.

        Args:
            genome_digest: The genome digest.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            seek_key_name: The name of the seek key.
            force_exists: Whether to raise an error if the path doesn't exist.
            abs_path: Return the content (``data/``) path rather than the alias path.

        Returns:
            str: The path to the asset file (for path-based seek keys), or the
                value directly (for non-path seek keys like string or json).
        """
        return self._seek(
            genome_digest,
            None,
            asset_group_name,
            asset_name,
            seek_key_name,
            force_exists=force_exists,
            abs_path=abs_path,
        )

    def seek_components(
        self,
        asset_registry_path_components: AssetRegistryPathComponents,
        force_exists: bool = False,
        abs_path: bool = False,
    ) -> str:
        """
        Get the path to an asset file from registry path components, or the value
        directly for non-path seek keys.

        The genome in a registry path is an alias. The default result is the
        alias-tree path built from exactly that alias and the asset name the
        caller supplied, with the seek value's genome-digest substring rewritten
        to the alias (matching the alias tree's per-file rename). ``abs_path``
        returns the digest-addressed content path instead.

        Args:
            asset_registry_path_components: The asset registry path components.
            force_exists: Whether to raise an error if the path doesn't exist.
            abs_path: Return the content (``data/``) path rather than the alias path.

        Returns:
            str: The path to the asset file (for path-based seek keys), or the
                value directly (for non-path seek keys like string or json).

        Raises:
            MissingAliasError: If the genome alias is not known here.
        """
        logger.debug(f"Seeking {asset_registry_path_components}")
        if (
            asset_registry_path_components.genome is None
            or asset_registry_path_components.asset_group is None
        ):
            raise ValueError("Genome name and asset name are required")
        genome_alias = GenomeAlias(asset_registry_path_components.genome)
        return self._seek(
            self._alias_manager.resolve(genome_alias),
            genome_alias,
            asset_registry_path_components.asset_group,
            asset_registry_path_components.asset,
            asset_registry_path_components.seek_key,
            force_exists=force_exists,
            abs_path=abs_path,
        )

    def _seek(
        self,
        genome_digest: GenomeDigest,
        genome_alias: GenomeAlias | None,
        asset_group_name: str,
        asset_name: str | None,
        seek_key_name: str | None,
        *,
        force_exists: bool,
        abs_path: bool,
    ) -> str:
        """Seek, keying the alias-tree path by ``genome_alias`` when one is given."""
        asset_name, seek_key = self.seek_key.resolve(
            genome_digest, asset_group_name, asset_name, seek_key_name
        )

        if not is_path_type(seek_key.type):
            # Non-path types: return the value directly as a string
            return seek_key.value

        if seek_key.asset.path is None:
            raise ValueError(
                f"Seek key {seek_key.name} has no path defined for asset {seek_key.asset}"
            )

        # Choose the alias tree key. When the caller named no alias, use a local
        # alias if the genome has one, otherwise return the content (``data/``)
        # path.
        alias_for_tree = genome_alias
        if alias_for_tree is None and not abs_path:
            local_aliases = self._alias_manager.get_for_genome(genome_digest)
            alias_for_tree = local_aliases[0] if local_aliases else None

        if abs_path or alias_for_tree is None:
            seek_key_path = (self._genome_folder / seek_key.asset.path / seek_key.value).absolute()
        else:
            # Alias-tree path, named for the alias the caller supplied. The alias
            # tree rewrites the genome digest to the alias in filenames, so the
            # seek value is rewritten to match.
            rewritten_value = seek_key.value.replace(genome_digest, alias_for_tree)
            seek_key_path = (
                self._alias_folder
                / alias_for_tree
                / asset_group_name
                / asset_name
                / rewritten_value
            ).absolute()

        if force_exists and not seek_key_path.exists():
            raise FileNotFoundError(f"Seek key path not found: {seek_key_path}")
        return str(seek_key_path)
