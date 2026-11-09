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

from refgenie.db.engine import Engine
from refgenie.db.queries import asset as queries
from refgenie.db.queries import asset_group as group_queries
from refgenie.db.tables import Asset, is_path_type
from refgenie.exceptions import MissingAssetError
from refgenie.logger import logger
from refgenie.models import AssetRegistryPathComponents, GenomeAlias, GenomeDigest, Page
from refgenie.managers.base import ResourceManager
from refgenie.managers.asset.alias_tree import AliasTree
from refgenie.managers.asset.content import AssetContentManager
from refgenie.managers.asset.group import AssetGroupManager
from refgenie.managers.asset.links import AssetLinkManager
from refgenie.managers.asset.seek_key import SeekKeyManager
from refgenie.managers.asset.tables import asset_table
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import Change
from refgenie.utils.search import Search

if TYPE_CHECKING:
    from refgenie.managers.alias import AliasManager
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.asset_class import AssetClassManager


class AssetManager(ResourceManager):
    """
    Manager for asset records: lookups, listing, removal, rename and local seek.

    Composed managers are public attributes: ``content`` (the write path: ``add``,
    ``adopt_name``), ``group`` (asset groups and their
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

    def get_by_digest(self, digest: str, *, with_links: bool = False) -> Asset:
        """
        Get an asset by its content digest.

        The asset comes with its group, genome, asset class, names and seek
        keys loaded; ``with_links`` also loads its parents and children, each
        loaded the same way.

        Args:
            digest: The digest of the asset.
            with_links: Also load the asset's parents and children.

        Returns:
            Asset: The asset.

        Raises:
            MissingAssetError: If no asset holds that digest.
        """
        with self._database_session as session:
            return queries.get_by_digest(session, digest, with_links=with_links)

    def exists_digest(self, digest: str) -> bool:
        """Whether some asset holds content ``digest``."""
        with self._database_session as session:
            return queries.find_by_digest(session, digest) is not None

    def page(
        self,
        *,
        genome_digest: GenomeDigest | None = None,
        asset_group_name: str | None = None,
        asset_group_id: int | None = None,
        name: str | None = None,
        recipe_name: str | None = None,
        search: Search | None = None,
        offset: int,
        limit: int,
    ) -> Page[Asset]:
        """
        One page of assets, by genome digest, group name and digest.

        Args:
            genome_digest: Only this genome's assets.
            asset_group_name: Only groups with this name.
            asset_group_id: Only this group.
            name: Only assets with this name (any of their names).
            recipe_name: Only assets built by this recipe.
            search: Match ``name``, ``digest`` or ``path``.
            offset: Rows to skip.
            limit: Most rows to return.

        Raises:
            ValueError: If the search names a field that is not searchable.
        """
        with self._database_session as session:
            rows, total = queries.list_page(
                session,
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_group_id=asset_group_id,
                name=name,
                recipe_name=recipe_name,
                search=search,
                offset=offset,
                limit=limit,
            )
        return Page(items=rows, offset=offset, limit=limit, total=total)

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
            return queries.get_by_name(session, genome_digest, asset_group_name, asset_name)

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

        Resolves the name to a digest and calls :meth:`remove_by_digest`.

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
            asset = queries.get_by_name(session, genome_digest, asset_group_name, asset_name)
            # Removal destroys every name for the content, so require the canonical one.
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
            asset = queries.get_by_digest(session, digest, with_links=True)

            if asset.children:
                raise ValueError(
                    f"Can't remove. Asset '{asset.registry_path}' has children. "
                    f"Remove them first: {', '.join(c.registry_path for c in asset.children)}"
                )

            genome_digest = asset.asset_group.genome.digest
            asset_group_name = asset.asset_group.name
            asset_names = [n.name for n in asset.asset_names]
            asset_digest = asset.digest
            registry_path = asset.registry_path

            canonical_name = asset.name
            session.delete(asset)
            session.commit()

            remaining = queries.list_all(session, [genome_digest], asset_group_name)

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

        # Resolve the group's directories first: removing the last group removes
        # the genome and the aliases these paths are keyed by.
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
        Check if an asset name exists. For content by digest, use :meth:`exists_digest`.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The genome digest.

        Returns:
            bool: Whether the asset exists. False for an unknown genome.
        """
        with self._database_session as session:
            return queries.exists_name(session, genome_digest, asset_group_name, asset_name)

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
        with self._database_session as session:
            return queries.list_all(session, genome_digests, asset_group_name)

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
            for name in sorted(n.name for n in asset.asset_names):
                asset_strings = (
                    [
                        f"{asset.asset_group.name}.{seek_key.name}:{name}"
                        for seek_key in asset.seek_keys
                    ]
                    if include_seek_keys
                    else [f"{asset.asset_group.name}:{name}"]
                )
                return_dict[asset.asset_group.genome.digest].extend(asset_strings)

        genome_digests_with_assets = list(return_dict.keys())
        aliases_dict = {}
        if genome_digests_with_assets:
            local_aliases = self._alias_manager.list_all()
            aliases_by_genome = defaultdict(list)
            for alias in local_aliases:
                if alias.genome_digest in genome_digests_with_assets:
                    aliases_by_genome[alias.genome_digest].append(alias.name)
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
            row = group_queries.find_name(session, asset_group.id, asset_name)
            if row is None:
                raise MissingAssetError(
                    genome=genome_digest, asset_group=asset_group_name, asset=asset_name
                )
            if group_queries.find_name(session, asset_group.id, new_asset_name) is not None:
                raise ValueError(
                    f"Asset name '{genome_digest}/{asset_group_name}:{new_asset_name}' "
                    f"already exists."
                )
            asset = queries.get_by_digest(session, row.asset_digest)
            # Content does not move. A canonical name also renames the asset, so
            # the staged-archive key follows.
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

        # Clear both trees the old name owned: a stale build directory holds the
        # completion flag, so a later build of the old name would be skipped.
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
            return seek_key.value

        # With no alias given, use the genome's first local alias, else the data/ path.
        alias_for_tree = genome_alias
        if alias_for_tree is None and not abs_path:
            local_aliases = self._alias_manager.get_for_genome(genome_digest)
            alias_for_tree = local_aliases[0] if local_aliases else None

        if abs_path or alias_for_tree is None:
            seek_key_path = (self._genome_folder / seek_key.asset.path / seek_key.value).absolute()
        else:
            # The alias tree renames the genome digest to the alias in filenames.
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
