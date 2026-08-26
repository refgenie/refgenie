"""
AssetManager - handles asset CRUD operations, path resolution, and defaults.

This is the single public interface for all asset operations. Internal classes
(AssetBuilder, AssetPuller, AssetRelations) are owned and managed by this class.
"""

import shutil
from pathlib import Path
from collections.abc import Callable, Iterable
from typing import Any, TYPE_CHECKING
from collections import defaultdict

from rich.table import Table
from sqlalchemy.engine import Engine
from sqlalchemy.orm import selectinload
from sqlmodel import select

from refgenie.db.tables import (
    Asset,
    AssetGroup,
    AssetName,
    Genome,
    is_path_type,
)
from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingAssetGroupError,
)
from refgenie.logger import logger
from refgenie.models import AssetRegistryPathComponents, BuildParams
from refgenie.managers.base import ResourceManager
from refgenie.managers.asset.alias_tree import AssetAliasTreeMixin
from refgenie.managers.asset.content import AssetContentMixin
from refgenie.managers.asset.queries import asset_by_name_stmt
from refgenie.managers.queries import one_or_raise
from refgenie.utils.tables import SECTION, build_table
from refgenie.managers.asset.seek import AssetSeekMixin
from refgenie.utils.prompt import Confirmer

if TYPE_CHECKING:
    from refgenie.managers.alias import AliasManager
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.asset_class import AssetClassManager
    from refgenie.managers.recipe import RecipeManager
    from refgenie.managers.sources.manager import SourceManager
    from refgenie.managers.asset.builder import AssetBuilder
    from refgenie.managers.asset.puller import AssetPuller
    from refgenie.managers.asset.relations import AssetRelations


class AssetManager(AssetContentMixin, AssetAliasTreeMixin, AssetSeekMixin, ResourceManager):
    """
    Manager for asset operations.

    This is the single public interface for all asset operations. Handles:
    - CRUD operations for assets
    - Build operations (via internal AssetBuilder)
    - Pull operations (via internal AssetPuller)
    - Parent/child relationships (via internal AssetRelations)
    """

    def __init__(
        self,
        database_engine: Engine,
        genome_folder: Path,
        alias_folder: Path,
        alias_manager: "AliasManager",
        genome_manager: "GenomeManager",
        asset_class_manager: "AssetClassManager",
        recipe_manager: "RecipeManager",
        source_manager: "SourceManager",
    ):
        """
        Initialize the AssetManager.

        Args:
            database_engine: The database engine.
            genome_folder: Path to genome data folder.
            alias_folder: Path to alias folder.
            alias_manager: The AliasManager for resolving genome names.
            genome_manager: The GenomeManager for genome operations.
            asset_class_manager: The AssetClassManager.
            recipe_manager: The RecipeManager for getting recipes.
            source_manager: The SourceManager for server clients.
        """
        super().__init__(database_engine)
        self._genome_folder = genome_folder
        self._alias_folder = alias_folder
        self._alias_manager = alias_manager
        self._genome_manager = genome_manager
        self._asset_class_manager = asset_class_manager
        self._recipe_manager = recipe_manager
        self._source_manager = source_manager
        # Internal instances (lazy init)
        self._builder: "AssetBuilder" | None = None
        self._puller: "AssetPuller" | None = None
        self._relations: "AssetRelations" | None = None

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

    # === Internal Lazy Properties ===

    @property
    def _asset_builder(self) -> "AssetBuilder":
        """Lazily initialize the internal AssetBuilder."""
        if self._builder is None:
            from refgenie.managers.asset.builder import AssetBuilder

            self._builder = AssetBuilder(
                database_engine=self.database_engine,
                genome_folder=self._genome_folder,
                alias_folder=self._alias_folder,
                data_folder=self.data_folder,
                recipe_manager=self._recipe_manager,
                asset_class_manager=self._asset_class_manager,
                asset_manager=self,
                alias_manager=self._alias_manager,
                asset_relations=self._asset_relations,
                genome_manager=self._genome_manager,
            )
        return self._builder

    @property
    def _asset_puller(self) -> "AssetPuller":
        """Lazily initialize the internal AssetPuller."""
        if self._puller is None:
            from refgenie.managers.asset.puller import AssetPuller

            self._puller = AssetPuller(
                database_engine=self.database_engine,
                genome_folder=self._genome_folder,
                alias_folder=self._alias_folder,
                source_manager=self._source_manager,
                asset_manager=self,
                alias_manager=self._alias_manager,
                genome_manager=self._genome_manager,
                asset_relations=self._asset_relations,
            )
        return self._puller

    @property
    def _asset_relations(self) -> "AssetRelations":
        """Lazily initialize the internal AssetRelations."""
        if self._relations is None:
            from refgenie.managers.asset.relations import AssetRelations

            self._relations = AssetRelations(self.database_engine)
        return self._relations

    def _resolve_genome_digest(self, *, genome_digest: str | None, genome_name: str | None) -> str:
        """
        Resolve genome digest from digest or name.

        Args:
            genome_digest: The digest of the genome.
            genome_name: The name of the genome.

        Returns:
            str: The resolved genome digest.
        """
        if genome_digest is None:
            if genome_name is None:
                raise ValueError("Either genome_digest or genome_name must be provided")
            genome_digest = self._resolve_alias_or_digest(genome_name)
        return genome_digest

    def _resolve_alias_or_digest(self, token: str) -> str:
        """
        Resolve a genome token that may be either an alias or a genome digest.

        Aliases take precedence: if ``token`` resolves as an alias, that digest
        wins. Otherwise, if ``token`` is itself a known genome digest, it is
        returned as-is. An unknown token re-raises ``MissingAliasError``.

        Args:
            token: An alias name or a genome digest.

        Returns:
            str: The resolved genome digest.

        Raises:
            MissingAliasError: If the token is neither a known alias nor a
                known genome digest.
        """
        try:
            return self._alias_manager.resolve(token)
        except MissingAliasError:
            if self._genome_manager.exists(token):
                return token
            raise

    def _resolve_genome_digests(
        self,
        genome_names: list[str] | None = None,
        genome_digests: list[str] | None = None,
        allow_missing: bool = False,
    ) -> list[str] | None:
        """
        Resolve genome digests from genome names or digests.

        Args:
            genome_names: The names of the genomes.
            genome_digests: The digests of the genomes.
            allow_missing: Whether to allow missing genomes.

        Returns:
            list[str]: The resolved genome digests.
        """
        if genome_digests is None:
            if genome_names is None:
                if not allow_missing:
                    raise ValueError("Either genome_digests or genome_names must be provided")
                return None
            genome_digests = [
                self._resolve_alias_or_digest(genome_name) for genome_name in genome_names
            ]
        return genome_digests

    # === CRUD Operations ===

    def get_by_digest(self, digest: str) -> Asset | None:
        """
        Get an asset by its content digest, or None if no asset holds that digest.

        Unlike ``get(digest=...)`` this does not raise when there is no match, so it
        can be used to test whether a piece of content is already in the catalog.

        Args:
            digest: The digest of the asset.

        Returns:
            Asset | None: The asset, or None if not found.
        """
        statement = (
            select(Asset)
            .where(Asset.digest == digest)
            .options(
                selectinload(Asset.seek_keys),
                selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
            )
        )
        with self._database_session as session:
            return session.exec(statement).unique().one_or_none()

    def get(
        self,
        asset_group_name: str | None = None,
        asset_name: str | None = None,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
        digest: str | None = None,
    ) -> Asset:
        """
        Get an asset by its name or digest.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.
            digest: The digest of the asset.

        Returns:
            Asset: The asset.
        """
        assert (asset_group_name is not None and asset_name is not None) or digest is not None, (
            "Either asset_group_name and asset_name or digest must be provided"
        )

        if digest is not None:
            statement = (
                select(Asset)
                .where(Asset.digest == digest)
                .options(
                    selectinload(Asset.seek_keys),
                    selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
                )
            )
            with self._database_session as session:
                return one_or_raise(
                    session, statement, MissingAssetError(digest=digest), unique=True
                )
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        statement = asset_by_name_stmt(
            genome_digest,
            asset_group_name,
            asset_name,
            options=[
                selectinload(Asset.seek_keys),
                selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
                selectinload(Asset.asset_names),
            ],
        )
        with self._database_session as session:
            return one_or_raise(
                session,
                statement,
                MissingAssetError(
                    genome=genome_digest, asset_group=asset_group_name, asset=asset_name
                ),
                unique=True,
            )

    def get_asset_dir(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ) -> Path:
        """
        Absolute path to the asset's data directory.

        This is the single definition of where an asset's files live. Callers
        that need the directory must use this rather than deriving it from a
        seek key, whose value may be nested and whose type may not be a path.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.

        Returns:
            Path: The absolute asset directory.
        """
        asset = self.get(
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            genome_name=genome_name,
            genome_digest=genome_digest,
        )
        if asset.path is None:
            raise ValueError(f"Incomplete asset, path is not set: {asset.registry_path}")
        return self._genome_folder / asset.path

    def get_group(
        self,
        asset_group_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ) -> AssetGroup:
        """
        Get an asset group by its name.

        Args:
            asset_group_name: The name of the asset group.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.

        Returns:
            AssetGroup: The asset group.
        """
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        statement = (
            select(AssetGroup)
            .join(Genome)
            .where(AssetGroup.name == asset_group_name, Genome.digest == genome_digest)
        )
        with self._database_session as session:
            return one_or_raise(
                session,
                statement,
                MissingAssetGroupError(genome=genome_digest, asset_group=asset_group_name),
                unique=True,
            )

    def remove(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
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
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.
            keep_asset_group: Whether to keep the asset group if no assets left.

        Returns:
            The registry path of the removed asset.

        Raises:
            ValueError: If ``asset_name`` is not the content's canonical name.
            MissingAssetError: If the asset does not exist.
        """
        genome_digest = self._resolve_genome_digest(
            genome_digest=genome_digest, genome_name=genome_name
        )

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

            session.delete(asset)
            session.commit()

            remaining = [a for a in asset_group.assets if a.digest != asset_digest]

        logger.info(f"Removed asset '{registry_path}'")
        self._purge_name_trees(genome_digest, asset_group_name, asset_names)

        if remaining or keep_asset_group:
            return registry_path

        # Empty group cleanup. The group's directories are resolved *before* the
        # group row goes away: removing the last group of a genome removes the
        # genome, and with it the aliases these paths are keyed by.
        logger.info(f"No assets left for '{genome_digest}/{asset_group_name}'. Removing group.")
        group_trees = self._owned_trees(genome_digest, asset_group_name)
        self.remove_group(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
        )
        for alias_name, path in group_trees:
            logger.info(f"Removing group files for alias '{alias_name}': {path}")
            shutil.rmtree(path, ignore_errors=True)

        return registry_path

    def remove_group(
        self,
        asset_group_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ):
        """
        Remove an asset group and all assets.

        Args:
            asset_group_name: The name of the asset group.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.
        """
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        with self._database_session as session:
            asset_group = self.get_group(
                genome_digest=genome_digest, asset_group_name=asset_group_name
            )
            session.delete(asset_group)
            session.commit()
            # check if there are any assets left for the genome and remove the genome if not
            if (
                not session.exec(
                    select(AssetGroup).join(Genome).where(Genome.digest == genome_digest)
                )
                .unique()
                .all()
            ):
                logger.info(
                    f"No assets groups left for the genome '{genome_digest}'. Removing the genome"
                )
                self._genome_manager.remove(genome_digest)
        logger.info(f"Removed asset group and all assets '{genome_digest}/{asset_group_name}'")

    def exists(
        self,
        asset_group_name: str | None = None,
        asset_name: str | None = None,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
        digest: str | None = None,
    ) -> bool:
        """
        Check if an asset exists.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.
            digest: The digest of the asset.

        Returns:
            bool: Whether the asset exists.
        """
        assert (asset_group_name is not None and asset_name is not None) or digest is not None, (
            "Either asset_group_name and asset_name or digest must be provided"
        )

        if digest is not None:
            statement = select(Asset).where(Asset.digest == digest)
            with self._database_session as session:
                result = session.exec(statement)
                return bool(result.first())

        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        statement = (
            select(AssetName)
            .join(AssetGroup, AssetGroup.id == AssetName.asset_group_id)
            .join(Genome)
            .where(
                AssetName.name == asset_name,
                AssetGroup.name == asset_group_name,
                Genome.digest == genome_digest,
            )
        )
        with self._database_session as session:
            result = session.exec(statement)
            return bool(result.first())

    def group_exists(
        self,
        asset_group_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ) -> bool:
        """
        Check if an asset group exists.

        Args:
            asset_group_name: The name of the asset group.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.

        Returns:
            bool: Whether the asset group exists.
        """
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        statement = (
            select(AssetGroup)
            .join(Genome)
            .where(AssetGroup.name == asset_group_name, Genome.digest == genome_digest)
        )
        with self._database_session as session:
            result = session.exec(statement)
            return bool(result.first())

    def list_assets(
        self,
        genome_digests: list[str] | None = None,
        genome_names: list[str] | None = None,
        asset_group_name: str | None = None,
    ) -> Iterable[Asset]:
        """
        List all assets.

        Args:
            genome_digests: The digests of genomes to filter by.
            genome_names: The names of genomes to filter by.
            asset_group_name: The name of asset group to filter by.

        Returns:
            Iterable[Asset]: A list of all assets.
        """
        genome_digests = self._resolve_genome_digests(
            genome_names=genome_names, genome_digests=genome_digests, allow_missing=True
        )
        statement = (
            select(Asset)
            .join(AssetGroup)
            .join(Genome)
            .options(selectinload(Asset.asset_group).selectinload(AssetGroup.genome))
            .options(selectinload(Asset.seek_keys))
            .options(selectinload(Asset.asset_names))
        )
        if genome_digests is not None:
            statement = statement.where(Genome.digest.in_(genome_digests))

        if asset_group_name is not None:
            statement = statement.where(AssetGroup.name == asset_group_name)

        with self._database_session as session:
            result = session.exec(statement)
            return result.unique().all()

    def list_seek_keys_values(
        self,
        genome_names: list[str] | None = None,
        genome_digests: list[str] | None = None,
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
            genome_names: Optional genome names to filter by.
            genome_digests: Optional genome digests to filter by.
            asset_group_name: Optional asset group name to filter by.

        Returns:
            A nested dict keyed by primary genome alias, asset group, asset,
            and seek key name with str values.
        """
        assets = self.list_assets(
            genome_names=genome_names,
            genome_digests=genome_digests,
            asset_group_name=asset_group_name,
        )

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

    def list_groups(
        self,
        genome_digests: list[str] | None = None,
        genome_names: list[str] | None = None,
    ) -> Iterable[AssetGroup]:
        """
        List all asset groups.

        Args:
            genome_digests: The digests of genomes to filter by.
            genome_names: The names of genomes to filter by.

        Returns:
            Iterable[AssetGroup]: A list of all asset groups.
        """
        genome_digests = self._resolve_genome_digests(
            genome_names=genome_names, genome_digests=genome_digests, allow_missing=True
        )
        statement = (
            select(AssetGroup)
            .join(Genome)
            .options(selectinload(AssetGroup.genome))
            .options(selectinload(AssetGroup.assets))
        )
        if genome_digests is not None:
            statement = statement.where(Genome.digest.in_(genome_digests))

        with self._database_session as session:
            result = session.exec(statement)
            return result.unique().all()

    def list_all(
        self,
        genome_names: list[str] | None = None,
        genome_digests: list[str] | None = None,
        asset_group_name: str | None = None,
        include_seek_keys: bool = False,
    ) -> tuple[dict[str, list[str]], dict[str, str]]:
        """
        List assets with formatted output.

        Args:
            genome_names: The names of genomes.
            genome_digests: The digests of genomes.
            asset_group_name: The name of the asset group.
            include_seek_keys: Whether to include seek keys.

        Returns:
            tuple[dict[str, list[str]], dict[str, str]]: Asset data and aliases data.
        """
        assets = self.list_assets(
            genome_digests=genome_digests,
            genome_names=genome_names,
            asset_group_name=asset_group_name,
        )
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
        genome_names: list[str] | None = None,
        genome_digests: list[str] | None = None,
        include_seek_keys: bool = False,
    ) -> list[Table]:
        """
        Create a table of all assets.

        Args:
            genome_names: The names of genomes.
            genome_digests: The digests of genomes.
            include_seek_keys: Whether to include seek keys.

        Returns:
            list[Table]: A list of tables.
        """
        asset_data, aliases_data = self.list_all(
            genome_names=genome_names,
            genome_digests=genome_digests,
            include_seek_keys=include_seek_keys,
        )
        asset_data_by_source = {"local": asset_data}
        aliases_data_by_source = {"local": aliases_data}

        return [
            self._create_table(
                asset_data=asset_data,
                aliases_data=aliases_data_by_source.get(source, {}),
                include_seek_keys=include_seek_keys,
                source=source,
            )
            for source, asset_data in asset_data_by_source.items()
        ]

    def remote_table(
        self,
        genome_digests: list[str] | None = None,
        include_seek_keys: bool = False,
        server_urls: list[str] | None = None,
    ) -> list[Table]:
        """
        Create a table of assets available on remote servers, one table per server.

        Args:
            genome_digests: The digests of genomes to filter by.
            include_seek_keys: Whether to include seek keys.
            server_urls: Optional list of server URLs to query. If not provided,
                         uses all subscribed servers.

        Returns:
            list[Table]: A list of tables, one per remote server.
        """
        asset_data_by_source, aliases_data_by_source = self.list_remote(
            genome_digests=genome_digests,
            include_seek_keys=include_seek_keys,
            server_urls=server_urls,
        )
        return [
            self._create_table(
                asset_data=asset_data,
                aliases_data=aliases_data_by_source.get(source, {}),
                include_seek_keys=include_seek_keys,
                source=source,
            )
            for source, asset_data in asset_data_by_source.items()
        ]

    def _create_table(
        self,
        asset_data: dict[str, list[str]],
        source: str,
        aliases_data: dict[str, str | None] = None,
        include_seek_keys: bool = False,
    ) -> Table:
        """
        Create a rich.Table from asset data.

        Args:
            asset_data: Asset data dictionary.
            source: The source of the data.
            aliases_data: Aliases data dictionary.
            include_seek_keys: Whether to include seek keys.

        Returns:
            Table: A Rich table.
        """
        if aliases_data is None:
            aliases_data = {}

        title = f"Refgenie assets. Source: {source}"
        if not asset_data:
            logger.warning("No assets found")
            return build_table(title, [], [])

        columns = ["Aliases", "Genome digest", "Asset group", "Asset"]
        if include_seek_keys:
            columns.append("Seek key")

        rows = []
        for genome_digest, asset_strings in asset_data.items():
            # One section per genome, so a long listing stays readable.
            if rows:
                rows.append(SECTION)
            aliases_str = aliases_data.get(genome_digest, "")
            for asset_string in asset_strings:
                components = AssetRegistryPathComponents.parse_registry_path(asset_string)
                row = [aliases_str, genome_digest, components.asset_group, components.asset]
                if include_seek_keys:
                    row.append(components.seek_key)
                rows.append(row)
        return build_table(title, columns, rows)

    def rename(
        self,
        asset_group_name: str,
        asset_name: str,
        new_asset_name: str,
        *,
        genome_digest: str | None = None,
        genome_name: str | None = None,
    ):
        """
        Rename an asset.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            new_asset_name: The new name for the asset.
            genome_digest: The digest of the genome.
            genome_name: The name of the genome.

        Returns:
            The renamed asset object.
        """
        genome_digest = self._resolve_genome_digest(
            genome_name=genome_name, genome_digest=genome_digest
        )
        asset_group = self.get_group(genome_digest=genome_digest, asset_group_name=asset_group_name)

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

        # Render the alias tree under the new name and clear the stale trees the
        # old name owned. Both of them: the build directory holds the completion
        # flag, and leaving it behind means a later build of the OLD name is
        # skipped as already done, while a rename back to it fails with "already
        # exists". _owned_trees is the same enumeration _purge_name_trees uses,
        # so the two removal paths cannot disagree about what a name owns.
        self._render_asset_alias_symlinks(genome_digest, asset_group_name, new_asset_name)
        for alias_name, path in self._owned_trees(genome_digest, asset_group_name, asset_name):
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
    # Delegations to internal collaborators: _asset_builder (build),
    # _asset_puller (pull / remote listing), _asset_relations (provenance)
    # =========================================================================

    def build(
        self,
        recipe_name: str,
        *,
        genome_name: str,
        genome_digest: str,
        asset_group_name: str,
        asset_name: str | None = None,
        recipe_version: str | None = None,
        params: BuildParams | None = None,
        docker: bool = False,
        docker_volumes: list[str] | None = None,
        asset_description: str | None = None,
        pipeline_kwargs: dict[str, Any | None] = None,
    ) -> Asset | None:
        """Build an asset using a recipe; returns the existing asset if it already exists."""
        return self._asset_builder.build(
            recipe_name=recipe_name, genome_name=genome_name, genome_digest=genome_digest,
            asset_group_name=asset_group_name, asset_name=asset_name,
            recipe_version=recipe_version, params=params, docker=docker,
            docker_volumes=docker_volumes, asset_description=asset_description,
            pipeline_kwargs=pipeline_kwargs,
        )

    def pull(
        self,
        asset_group_name: str,
        alias_name: str | None = None,
        genome_digest: str | None = None,
        asset_name: str | None = None,
        force: bool | None = None,
        force_large: bool | None = None,
        force_server_urls: list[str] | None = None,
        size_cutoff: int | float | None = None,
        sigint_handler: Callable | None = None,
        confirm: Confirmer | None = None,
    ) -> Asset | None:
        """Pull an asset from remote servers; returns the added asset, or None if pull failed."""
        return self._asset_puller.pull(
            asset_group_name=asset_group_name, alias_name=alias_name,
            genome_digest=genome_digest, asset_name=asset_name, force=force,
            force_large=force_large, force_server_urls=force_server_urls,
            size_cutoff=size_cutoff, sigint_handler=sigint_handler, confirm=confirm,
        )

    def list_remote(
        self,
        genome_digests: list[str] | None = None,
        include_seek_keys: bool = False,
        server_urls: list[str] | None = None,
    ) -> tuple[dict[str, dict[str, list[str]]], dict[str, dict[str, str]]]:
        """List assets on remote servers as (assets by server URL, aliases by server URL)."""
        return self._asset_puller.list_remote(
            genome_digests=genome_digests,
            include_seek_keys=include_seek_keys,
            server_urls=server_urls,
        )

    def list_remote_assets_for_genome(
        self,
        genome_digest: str,
        server_urls: list[str] | None = None,
    ) -> list[dict]:
        """List all assets available for a genome on remote servers."""
        return self._asset_puller.list_remote_assets_for_genome(
            genome_digest=genome_digest, server_urls=server_urls
        )

    def list_remote_genomes(
        self,
        server_urls: list[str] | None = None,
    ) -> list[dict]:
        """List all genomes available on remote servers."""
        return self._asset_puller.list_remote_genomes(server_urls=server_urls)

    def estimate_pull_size(
        self,
        asset_list: list[dict],
    ) -> tuple[int, int]:
        """Calculate (total_bytes, asset_count) for a list of remote asset dicts."""
        return self._asset_puller.estimate_pull_size(asset_list)

    def pull_multiple(
        self,
        asset_list: list[dict],
        force: bool | None = None,
        force_large: bool | None = None,
        size_cutoff: int | float | None = None,
        sigint_handler: Callable | None = None,
        confirm: Confirmer | None = None,
    ) -> list:
        """Pull multiple assets (pull() for each); returns the successfully pulled Assets."""
        return self._asset_puller.pull_multiple(
            asset_list=asset_list, force=force, force_large=force_large,
            size_cutoff=size_cutoff, sigint_handler=sigint_handler, confirm=confirm,
        )

    def init_genome_from_remote(
        self,
        alias_name: str,
        genome_digest: str | None = None,
        genome_description: str | None = None,
        server_urls: list[str] | None = None,
    ) -> bool:
        """Register a genome locally from remote metadata; False if not found or already exists."""
        return self._asset_puller.init_genome_from_remote(
            alias_name=alias_name, genome_digest=genome_digest,
            genome_description=genome_description, server_urls=server_urls,
        )

    def get_parents(
        self, genome_digest: str, asset_group_name: str, asset_name: str
    ) -> list[Asset]:
        """Get parent assets."""
        return self._asset_relations.get_parents(genome_digest, asset_group_name, asset_name)

    def get_children(
        self, asset_group_name: str, asset_name: str, genome_digest: str
    ) -> list[Asset]:
        """Get child assets."""
        return self._asset_relations.get_children(asset_group_name, asset_name, genome_digest)

    def get_size(self, asset_group_name: str, asset_name: str, genome_digest: str) -> int:
        """Get asset size in bytes (sum of its seek key sizes)."""
        return self._asset_relations.get_size(asset_group_name, asset_name, genome_digest)

    def set_parents(
        self,
        genome_digest: str,
        asset_group_name: str,
        asset_name: str,
        parent_asset_digests: list[str],
    ) -> None:
        """Set parent assets."""
        self._asset_relations.set_parents(
            genome_digest, asset_group_name, asset_name, parent_asset_digests
        )
