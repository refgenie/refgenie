"""
Refgenie is a SQLModel backed reference genome manager.
"""

import os
from contextlib import contextmanager
from pathlib import Path
from collections.abc import Callable, Generator
from typing import Any

from refget.store import RefgetStore
from sqlalchemy.engine import Engine as SqlalchemyDatabaseEngine
from sqlmodel import Session

from refgenie.config import config
from refgenie.config.settings import LogLevel
from refgenie.const import ALIAS_DIR, DEFAULT_SERVER_URL
from refgenie.core.bulk import BulkTransferMixin
from refgenie.core.lifecycle import DatabaseLifecycleMixin
from refgenie.core.populate import RegistryPathPopulateMixin
from refgenie.core.sequences import SequenceAccessMixin
from refgenie.db.events import register_events
from refgenie.db.tables import (
    Alias,
    Asset,
    Recipe,
    SeekKeyType,
)
from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingAssetGroupError,
    MissingBuildInputError,
    MissingGenomeError,
    MissingRecipeError,
)
from refgenie.logger import logger
from refgenie.managers import (
    AliasManager,
    StageManager,
    AssetClassManager,
    AssetManager,
    ConfigurationManager,
    GenomeManager,
    RecipeManager,
)
from refgenie.managers.alias import FederatedAliasManager
from refgenie.managers.store import StoreManager
from refgenie.managers.sources import ServerClient, SourceManager, make_source
from refgenie.core.mode import LocalMode, RefgenieMode, ServerMode
from refgenie.core.store_router import RefgetStoreRouter
from refgenie.models import (
    AssetRegistryPathComponents,
    BuildCommandValues,
    BuildParams,
)
from refgenie.utils.build import get_build_dir
from refgenie.utils.templating import jinja_render_template_strictly
from refgenie.utils.prompt import Confirmer, resolve_confirmer
from refgenie.utils.symlinks import (
    create_alias_symlinks,
    get_build_paths as get_build_paths_util,
    get_symlink_paths as get_symlink_paths_util,
)


class Refgenie(
    DatabaseLifecycleMixin,
    BulkTransferMixin,
    SequenceAccessMixin,
    RegistryPathPopulateMixin,
):
    """
    A reference genome manager.

    Architecture note: Do NOT add wrapper methods that simply delegate to a single
    manager (e.g., `get_asset()` wrapping `self.asset.get()`). Instead, callers
    should use managers directly via properties: `rgc.alias`, `rgc.genome`,
    `rgc.asset`, etc. The Refgenie class should only contain methods that require
    cross-manager coordination or provide functionality beyond what a single
    manager offers.
    """

    def __init__(
        self,
        database_config_path: str | Path | None = None,
        database_engine: SqlalchemyDatabaseEngine | None = None,
        server_clients_mapping: dict[str, ServerClient | None] = None,
        suppress_migrations: bool = False,
        server_mode: bool = False,
    ):
        """
        Initialize the reference genome manager.

        Args:
            database_config_path: The path to the database configuration file. If not provided,
                the default configuration is used.
            database_engine: The database engine to use. If not provided, the default engine is used.
            server_clients_mapping: A mapping of server URLs to server clients. If not provided,
                the default server clients are used.
            suppress_migrations: Whether to suppress migrations. If True, migrations are not applied even if needed.
            server_mode: Run in server mode (federated remote stores from the
                ``store`` registry table, SQL aliases, no local sequence
                ingestion). Local mode (the default) owns one on-disk store.
        """
        register_events()  # Register SQLAlchemy event handlers
        logger.debug(f"{config=}")
        self._database_engine = database_engine or self.get_default_database_engine(
            database_config_path=database_config_path,
            echo=config.log_level == LogLevel.DEBUG,
        )
        if server_clients_mapping is not None:
            # validate server clients mapping. The clients need to follow the ServerClient Protocol
            for url, client in server_clients_mapping.items():
                if not isinstance(client, ServerClient):
                    raise ValueError(
                        f"Invalid server client for {url}. Does not match ServerClient Protocol"
                    )
        self._recipe_manager = RecipeManager(self.database_engine)
        self._asset_class_manager = AssetClassManager(self.database_engine)
        self._configuration_manager = ConfigurationManager(self.database_engine)
        self._stage_manager = StageManager(self.database_engine)
        self._source_manager = SourceManager(
            database_engine=self.database_engine,
            configuration_manager=self._configuration_manager,
            server_clients=server_clients_mapping,
        )
        # Mode is determined once, here, and never checked again
        if server_mode:
            self._mode: RefgenieMode = ServerMode(database_engine=self._database_engine)
        else:
            self._mode = LocalMode(
                genome_folder_getter=lambda: self.genome_folder,
                database_engine=self._database_engine,
            )

        self._store_manager = StoreManager(self.database_engine)
        self._alias_manager = self._mode.create_alias_manager()
        self._sequences_enabled = self._mode.sequences_enabled
        self._store_router: RefgetStoreRouter | None = None  # Lazy, built by mode

        # Only local mode reads a store, and the store is built lazily, so the
        # getter is handed over here rather than at construction.
        if isinstance(self._alias_manager, FederatedAliasManager):
            self._alias_manager.set_store_getter(lambda: self.refget_store)
        self._genome_manager = GenomeManager(
            self.database_engine,
            refget_store_getter=lambda: self.refget_store,
            alias_manager_getter=lambda: self._alias_manager,
        )
        self._asset_manager: AssetManager | None = None  # Lazy init after migrations

        if not suppress_migrations and self.check_for_db_migrations():
            self.migrate_db()

    def __str__(self):
        return f"Refgenie(database_engine={self.database_engine}, server_clients={list(self.server_clients.values())})"

    def __repr__(self):
        return str(self)

    @property
    def store_router(self) -> RefgetStoreRouter:
        """The federated store router.

        Built lazily (after migrations) from the mode: one on-disk store in
        local mode, or every enabled ``Store`` row in server mode. Each backend
        loads collection/alias metadata only -- sequence bytes are fetched on
        demand in getseq().
        """
        if self._store_router is None:
            self._store_router = self._mode.create_store_router(self._store_manager)
        return self._store_router

    def reload_store_router(self) -> RefgetStoreRouter:
        """Rebuild the router from the current ``store`` registry.

        Call after ``store add``/``sync``/``remove`` so a long-lived instance
        picks up registry changes without a restart. The alias manager caches
        against the store it reads, so it is invalidated with the router.
        """
        self._store_router = None
        self._alias_manager.invalidate()
        return self.store_router

    @property
    def refget_store(self) -> RefgetStore:
        """The default (highest-priority / writable) store in the router.

        Write paths (local-mode genome init, FHR sidecars, store-backed aliases)
        operate on this single store. Genome-specific read paths should route
        through :attr:`store_router` instead so federation dispatches correctly.
        """
        return self.store_router.default_store

    @property
    @contextmanager
    def _database_session(self) -> Generator[Session, None, None]:
        """
        Provide a transactional scope around a series of query
        operations.
        """
        session = Session(self.database_engine, expire_on_commit=False)
        try:
            yield session
        except:
            session.rollback()
            raise
        finally:
            session.close()

    @property
    def recipe(self) -> RecipeManager:
        """The recipe manager."""
        return self._recipe_manager

    @property
    def asset_class(self) -> AssetClassManager:
        """The asset class manager."""
        return self._asset_class_manager

    @property
    def configuration(self) -> ConfigurationManager:
        """The configuration manager."""
        return self._configuration_manager

    @property
    def stage(self) -> StageManager:
        """The stage manager."""
        return self._stage_manager

    @property
    def sources(self) -> SourceManager:
        """The source manager for external data sources."""
        return self._source_manager

    @property
    def store(self) -> StoreManager:
        """The store manager for the federation registry."""
        return self._store_manager

    @property
    def alias(self) -> AliasManager:
        """The alias manager."""
        return self._alias_manager

    @property
    def genome(self) -> GenomeManager:
        """The genome manager."""
        return self._genome_manager

    @property
    def asset(self) -> AssetManager:
        """
        Get the asset manager.

        This is the single public interface for all asset operations.

        Returns:
            AssetManager: The asset manager.
        """
        if self._asset_manager is None:
            self._asset_manager = AssetManager(
                database_engine=self.database_engine,
                genome_folder=self.genome_folder,
                alias_folder=self.alias_folder,
                alias_manager=self._alias_manager,
                genome_manager=self._genome_manager,
                asset_class_manager=self._asset_class_manager,
                recipe_manager=self._recipe_manager,
                source_manager=self._source_manager,
            )
        return self._asset_manager

    @property
    def database_engine(self) -> SqlalchemyDatabaseEngine:
        """
        Get the database engine.

        Returns:
            Engine: The database engine.
        """
        return self._database_engine

    @property
    def server_clients(self) -> dict[str, ServerClient]:
        """
        Get the server clients mapping.

        Returns:
            dict[str, ServerClient]: The server clients mapping.
        """
        return self._source_manager.server_clients

    @property
    def genome_folder(self) -> Path:
        """
        Get the genomes directory.

        Returns:
            Path: The genomes directory.
        """
        configuration = self.configuration.get_latest()
        return Path(configuration.genome_folder)

    @property
    def genome_stage_folder(self) -> Path | None:
        """
        Get the genome stage directory, if not set, return None.

        Returns:
            Path: The genome stage directory.
        """
        configuration = self.configuration.get_latest()
        return None if (a := configuration.genome_stage_folder) is None else Path(a)

    @property
    def data_folder(self) -> Path:
        """
        Get the data directory.

        Returns:
            Path: The data directory.
        """
        return self.genome_folder / "data"

    @property
    def alias_folder(self) -> Path:
        """
        Get the alias directory.

        Returns:
            Path: The alias directory.
        """
        return self.genome_folder / ALIAS_DIR

    def get_genome_init_target_template(self) -> Path:
        """Return a path template for genome init sentinel files.

        The returned path contains ``{genome_name}`` for Snakemake wildcard
        substitution.  After ``refgenie1 genome init`` succeeds the sentinel
        file is touched so that downstream Snakemake rules can depend on it.
        """
        return self.alias_folder / "{genome_name}" / ".genome_init_complete"

    def get_asset_build_target_template(self, asset_group_name: str, asset_name: str) -> Path:
        """Return a Snakemake target path template for a built asset.

        The returned path contains ``{genome_name}`` for Snakemake wildcard
        substitution. It is the build-stats flag file that ``refgenie build``
        writes on success, used by the generated Snakefile to wire the DAG.
        Public wrapper around the asset builder's internal target template.
        """
        return self.asset._asset_builder._get_build_target_template(
            asset_group_name=asset_group_name, asset_name=asset_name
        )

    @classmethod
    def parse_asset_registry_path(cls, asset_registry_path: str) -> AssetRegistryPathComponents:
        """
        Split an asset registry path into its components.

        For example 'genome_name/asset_group_name.seek_key_name:asset_name'
        to: ('genome_name', 'asset_group_name', 'seek_key_name', 'asset_name')

        Args:
            asset_registry_path: The asset registry path.

        Returns:
            AssetRegistryPathComponents: The genome, asset group, seek key, and asset names.
        """
        return AssetRegistryPathComponents.parse_registry_path(asset_registry_path)

    def add(
        self,
        asset_class_name: str,
        path: Path,
        asset_group_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
        asset_name: str | None = None,
        description: str | None = None,
        recipe: Recipe | None = None,
        custom_seek_keys: dict[str, str | tuple[str, "SeekKeyType"]] | None = None,
    ) -> Asset | None:
        """
        Register an existing asset from a path in the catalog.

        Args:
            genome_name: The name of the genome. Either genome_name or genome_digest must be provided.
            genome_digest: The digest of the genome. Either genome_name or genome_digest must be provided.
            asset_class_name: The name of the asset class.
            path: The path to the asset.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            description: The description of the asset.
            recipe: The recipe used to build the asset.
            custom_seek_keys: Dict mapping seek key names to values. Values can be
                a plain string (defaults to SeekKeyType.string) or a
                (value, SeekKeyType) tuple for explicit type control.

        Returns:
            Asset: The added asset.
        """
        asset = self.asset.add_from_path(
            asset_class_name=asset_class_name,
            path=path,
            asset_group_name=asset_group_name,
            genome_name=genome_name,
            genome_digest=genome_digest,
            asset_name=asset_name,
            description=description,
            recipe=recipe,
            custom_seek_keys=custom_seek_keys,
        )
        # Post-operation: create symlinks (like set_genome_alias does).
        # The genome may be identified by name OR by digest -- both are
        # documented as valid, so both must produce an alias-tree entry.
        if asset is not None:
            digest = self._resolve_genome_digest(
                genome_digest=genome_digest, genome_name=genome_name
            )
            aliases = self.alias.get_for_genome(digest)
            if aliases:
                # No whitelist: every alias of this genome gets its own tree.
                self._symlink_alias(
                    alias_name=aliases[0],
                    asset_group_name=asset_group_name,
                    asset_name=asset_name,
                )
            else:
                logger.warning(
                    f"Genome '{digest}' has no aliases; skipping alias tree for "
                    f"'{asset_group_name}:{asset_name}'"
                )
        return asset

    def _resolve_genome_digest(self, *, genome_digest: str | None, genome_name: str | None) -> str:
        """
        Most of the functions require genome digest to perform operations. This function
        resolves the genome digest from the provided digest or name.
        - If both are provided, the digest is used.
        - If only the name is provided, the digest is resolved from the name.
        - If only the digest is provided, the digest is used.

        Args:
            genome_digest: The digest of the genome.
            genome_name: The name of the genome.

        Returns:
            str: The resolved genome digest.
        """
        if genome_digest is None:
            if genome_name is None:
                raise ValueError("Either genome_digest or genome_name must be provided")
            genome_digest = self.alias.resolve(genome_name)
        return genome_digest

    def set_genome_alias(
        self,
        alias_name: str,
        genome_digest: str | None = None,
        genome_description: str | None = None,
        server_urls: list[str] | None = None,
    ) -> Alias:
        """
        Set a genome alias, possibly by querying the server for the digest and description.

        If the genome digest is not provided, the server is queried for the digest and description

        Args:
            alias_name: The name of the alias.
            genome_digest: The digest of the genome. Optional.
            genome_description: The description of the genome. Optional, and only used if the digest is provided.
            server_urls: The URLs of the server. Optional, and only used if the digest is not provided.

        Returns:
            Alias: The added alias.

        Raises:
            MissingAliasError: genome_digest not provided and the alias is not found on the server
        """
        # Step 1: make sure the genome exists locally. This method owns alias
        # registration (step 2), so never hand `alias_names` to the genome layer.
        if genome_digest is None:
            servers = server_urls or list(self.sources.get_subscriptions())
            source = None

            for server_url in servers:
                try:
                    source = make_source(server_url)
                except Exception:
                    continue
                resolved_digest = source.resolve_alias(alias_name)
                if resolved_digest is not None:
                    genome_digest = resolved_digest
                    logger.info(f"Resolved alias '{alias_name}' to digest: {genome_digest}")
                    break

            if genome_digest is None:
                raise MissingAliasError(alias_name)

            logger.info(f"Determined digest for {alias_name}: {genome_digest}")

            self.genome.initialize_genome(
                source=source,
                digest=genome_digest,
                description=genome_description or "",
                alias_names=[],
                use_existing=True,
            )
        else:
            try:
                self.genome.get(genome_digest)
            except MissingGenomeError:
                self.genome.add(
                    genome_digest,
                    genome_description or "No description provided",
                    [],
                )

        # Step 2: register the alias and render its view of the asset tree.
        alias = self._alias_manager.add(alias_name, genome_digest)
        self._symlink_alias(alias_name=alias_name, alias_whitelist=[alias_name])
        return alias

    def initialize_and_build(
        self,
        fasta_file_path: Path,
        genome_names: list[str],
        description: str = "",
        species_name: str | None = None,
        use_existing: bool = False,
        build_fasta: bool = True,
    ) -> tuple[str, bool]:
        """Initialize a genome from FASTA and optionally build the fasta asset.

        This is the recommended way to add a new genome. It ingests the FASTA
        into the RefgetStore and builds the fasta asset (fa, fai, chrom.sizes)
        in one step.

        Args:
            fasta_file_path: Path to the FASTA file.
            genome_names: Alias names for the genome.
            description: Genome description.
            species_name: Species name.
            use_existing: Allow re-init of existing genome.
            build_fasta: Whether to build the fasta asset (default True).

        Returns:
            Tuple of (genome_digest, was_created).
        """
        digest, created = self.genome.initialize_genome(
            fasta_file_path=fasta_file_path,
            description=description,
            alias_names=genome_names,
            use_existing=use_existing,
            species_name=species_name,
        )

        if build_fasta:
            self.build_asset(
                recipe_name="fasta",
                genome_name=genome_names[0],
                asset_group_name="fasta",
            )

        return digest, created

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
            genome_digests = [self.alias.resolve(genome_name) for genome_name in genome_names]
        return genome_digests

    @staticmethod
    def resolve_custom_seek_keys(custom_seek_keys: dict[str, str]) -> dict[str, Any]:
        """Resolve custom seek keys by executing shell commands.

        Used by the snakefile template to resolve recipe custom seek keys
        outside of a full build context.
        """
        from subprocess import check_output

        return (
            {
                key: check_output(cmd, shell=True).decode("utf-8").strip()
                for key, cmd in custom_seek_keys.items()
            }
            if custom_seek_keys
            else {}
        )

    @staticmethod
    def resolve_default_asset(default_asset: str, namespaces) -> str:
        """Resolve the default asset name using template rendering.

        Used by the snakefile template to determine the default asset name
        from a recipe's default_asset template string.
        """
        # Deliberately NO `or "default"` fallback. Assets are named after the
        # version of the tool that built them, and that version is produced by a
        # shell pipeline in the recipe's custom_seek_keys. Those pipelines end in
        # grep/awk/head, so a tool that fails or stops reporting its version
        # yields an empty string with a zero exit status -- silently.
        #
        # Never fall back to "default": it would bake a mis-named asset into
        # build paths and published keys.
        #
        # Recipes that legitimately have no tool version (fasta, fasta_index)
        # declare the literal `default_asset: "default"`, which renders truthy
        # and is unaffected.
        resolved = jinja_render_template_strictly(default_asset, namespaces)
        if not resolved or not resolved.strip():
            raise ValueError(
                f"default_asset template {default_asset!r} resolved to an empty asset "
                f"name. Refusing to fall back to 'default': that would build and "
                f"publish a mis-named asset. Check that the recipe's custom_seek_keys "
                f"commands run in this environment and actually print a version."
            )
        return resolved

    def build_asset(
        self,
        recipe_name: str,
        genome_name: str,
        asset_group_name: str,
        asset_name: str | None = None,
        recipe_version: str | None = None,
        params: BuildParams | None = None,
        stage: bool = False,
        docker: bool = False,
        docker_volumes: list[str] | None = None,
        asset_description: str | None = None,
        pull_parents: bool = False,
        pipeline_kwargs: dict[str, Any | None] = None,
        push_to: list[str] | None = None,
    ) -> Asset | None:
        """
        Build an asset for a specified genome based on provided recipe.

        This is a coordinator method that resolves the genome digest and
        delegates the actual build to AssetBuilder. The genome must already
        be initialized (via genome.initialize_genome or genome init CLI).

        Args:
            recipe_name: The name of the recipe to use
            genome_name: The name of the genome to build the asset for
            asset_group_name: The name of the asset group to use
            asset_name: The name of the asset to build
            recipe_version: The version of the recipe to use. If not provided, the latest version is used.
            params: The build parameters. Some recipes require input assets to be built.
            stage: Whether to stage the asset after building. Only possible if the stage folder is set.
            docker: Whether to use docker to build the asset.
            docker_volumes: The docker volumes to mount.
            asset_description: The description of the asset, used only if the asset does not exist.
            pull_parents: Whether to pull the parent genome if not found.
            pipeline_kwargs: Additional kwargs for PipelineManager.
            push_to: Optional list of remote names/IDs to create push intent records for.
        Returns:
            Asset: The asset that was built, or None if skipped/failed
        """
        # Fail before the build, not after it. A misconfigured stage folder used
        # to surface only once the build had finished, which from a browser (or
        # a multi-hour bowtie2 index) means losing the whole run to a check that
        # costs nothing up front.
        if stage and self.genome_stage_folder is None:
            raise ValueError(
                "Can't stage the asset: genome_stage_folder is not set. "
                "Set it in the refgenie config or build without --stage."
            )

        # Resolve genome digest, with optional pull_parents fallback
        try:
            genome_digest = self.alias.resolve(genome_name)
        except MissingAliasError:
            msg = f"Genome '{genome_name}' not found."
            if pull_parents:
                msg += " Attempting to pull from remote server."
                logger.warning(msg)
                self.pull(asset_group_name="fasta", alias_name=genome_name)
                genome_digest = self.alias.resolve(genome_name)
            else:
                msg += " Initialize it first with 'refgenie genome init'."
                logger.error(msg)
                raise

        # Delegate build to AssetManager
        asset = self.asset.build(
            recipe_name=recipe_name,
            genome_name=genome_name,
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            recipe_version=recipe_version,
            params=params,
            docker=docker,
            docker_volumes=docker_volumes,
            asset_description=asset_description,
            pipeline_kwargs=pipeline_kwargs,
        )

        # Stage if requested (Refgenie owns StageManager, not AssetManager)
        if stage and asset is not None:
            # genome_stage_folder was validated at the top of this method.
            self._stage_manager.create(
                asset=asset,
                genome_folder=self.genome_folder,
                genome_stage_folder=self.genome_stage_folder,
                push_to=push_to,
                build_dir=get_build_dir(
                    genome_folder=self.genome_folder,
                    genome_name=genome_name,
                    asset_group_name=asset_group_name,
                    asset_name=asset.name,
                ),
            )

        return asset

    def preflight_build(
        self,
        recipe_name: str,
        genome_name: str,
        asset_group_name: str,
        asset_name: str | None = None,
        recipe_version: str | None = None,
        params: BuildParams | None = None,
        stage: bool = False,
    ) -> dict[str, Any]:
        """Validate a prospective build without starting it.

        The seam the web build form's preflight endpoint calls: it runs the
        same resolution and validation `build_asset` would (genome, recipe,
        input assets, required files/params, staging config, default asset
        name) but starts no pipeline and writes nothing. Callers outside this
        class must use this method rather than reaching into
        ``AssetBuilder._``-private validators.

        Args:
            recipe_name: The name of the recipe to use.
            genome_name: The name (alias) of the genome to build for.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset to build. If omitted, the
                recipe's ``default_asset`` template is resolved.
            recipe_version: The version of the recipe. Latest if omitted.
            params: The build parameters (assets/params/files).
            stage: Whether the build would stage the asset afterwards.

        Returns:
            ``{"ok": bool, "errors": [{"field", "code", "message"}, ...],
            "resolved": {...}}`` -- ``errors`` is field-scoped so a form can
            attach each problem to its input; ``resolved`` reports what the
            build would use (genome digest, recipe version, input assets,
            asset name).
        """
        errors: list[dict[str, str]] = []
        resolved: dict[str, Any] = {}

        def problem(field: str, code: str, message: str) -> None:
            errors.append({"field": field, "code": code, "message": message})

        if stage and self.genome_stage_folder is None:
            problem(
                "stage",
                "missing_build_input",
                "Staging is not configured: genome_stage_folder is not set. "
                "Set it in the refgenie config or build without staging.",
            )

        genome_digest: str | None = None
        try:
            genome_digest = self.alias.resolve(genome_name)
            resolved["genome_digest"] = genome_digest
        except MissingAliasError:
            problem(
                "genome",
                "genome_not_found",
                f"Genome '{genome_name}' not found. Initialize it first with "
                "'refgenie genome init', or build with pull_parents.",
            )

        recipe = None
        try:
            recipe = self.recipe.get(recipe_name, recipe_version)
            resolved["recipe"] = recipe_name
            resolved["recipe_version"] = recipe.version
        except MissingRecipeError as exc:
            problem("recipe", "recipe_not_found", str(exc))

        if genome_digest is None or recipe is None:
            return {"ok": False, "errors": errors, "resolved": resolved}

        builder = self.asset._asset_builder
        if params is not None:
            params.populate_with_defaults_from_recipe(recipe)

        input_assets = None
        try:
            input_assets = builder._resolve_input_assets(
                genome_digest=genome_digest,
                recipe_name=recipe_name,
                recipe_version=recipe_version,
                build_params=params,
            )
            resolved["input_assets"] = {
                name: (asset.registry_path if asset is not None else None)
                for name, asset in (input_assets or {}).items()
            }
        except (MissingAssetError, MissingAssetGroupError, MissingAliasError) as exc:
            problem("params.assets", "asset_not_found", str(exc))
        except ValueError as exc:
            problem("params.assets", "conflict", str(exc))

        command_values = BuildCommandValues(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            params=params.params if params is not None else {},
            files=params.files if params is not None else {},
            assets=input_assets,
            custom_seek_keys={},
            genome_folder=self.genome_folder,
            refget_store_path=str(self.genome_folder / ".refget_store"),
        )
        try:
            builder._validate_build_inputs(recipe, command_values)
        except MissingBuildInputError as exc:
            problem("params", "missing_build_input", str(exc))

        if asset_name:
            resolved["asset_name"] = asset_name
        else:
            # The default-asset template needs the recipe's custom seek keys
            # (tool versions, produced by shell one-liners). Failing here is
            # itself useful preflight information: the tool is not installed.
            try:
                command_values.custom_seek_keys = (
                    builder._resolve_custom_seek_keys(recipe.custom_seek_keys or {}) or {}
                )
                resolved["asset_name"] = builder._resolve_default_asset(
                    recipe.default_asset, command_values
                )
            except Exception as exc:  # noqa: BLE001 - report, don't crash a preflight
                problem(
                    "asset",
                    "refgenie_error",
                    f"Could not resolve the default asset name: {exc}",
                )

        return {"ok": not errors, "errors": errors, "resolved": resolved}

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
        """
        Download and unpack an asset for a given reference genome.

        Args:
            asset_group_name: Name of a group of assets to fetch.
            alias_name: Name of a reference genome assembly of interest.
            genome_digest: Digest of the genome.
            asset_name: Name of particular asset to fetch.
            force: How to handle case in which asset path already exists.
            force_large: How to handle archives larger than size_cutoff (default 10GB).
            force_server_urls: Force specific server URLs to use.
            size_cutoff: Maximum archive file size to download without prompt.
            sigint_handler: Signal handler for interrupts during download.
            confirm: Confirmation callback. Defaults to a refusal unless the CLI
                has enabled interactive prompts; see `refgenie.utils.prompt`.

        Returns:
            Asset or None: The added asset, or None if pull failed.
        """
        # Check for subscriptions before pulling (subscribe prompt lives here, not in AssetPuller)
        if not force_server_urls and not self._source_manager.get_subscriptions():
            logger.error("No server subscriptions found")
            if not resolve_confirmer(confirm)("Would you like to subscribe to the default server?"):
                logger.info("Skipping pull")
                return None
            default_servers = [DEFAULT_SERVER_URL]
            self.configuration.subscribe(default_servers)

        # Delegate to AssetManager
        return self.asset.pull(
            asset_group_name=asset_group_name,
            alias_name=alias_name,
            genome_digest=genome_digest,
            asset_name=asset_name,
            force=force,
            force_large=force_large,
            force_server_urls=force_server_urls,
            size_cutoff=size_cutoff,
            sigint_handler=sigint_handler,
            confirm=confirm,
        )

    def get_symlink_paths(
        self,
        genome_digest: str,
        asset_group_name: str | None = None,
        asset_name: str | None = None,
        alias_whitelist: list[str] | None = None,
    ) -> dict[str, Path]:
        """
        Get alias-directory paths for the selected genome/asset group/asset,
        keyed by alias name.

        Args:
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            alias_whitelist: The list of aliases to include.

        Returns:
            dict[str, Path]: Mapping of alias names to their directory paths.
        """
        aliases = self._resolve_aliases(genome_digest, alias_whitelist)
        return get_symlink_paths_util(
            alias_folder=self.alias_folder,
            aliases=aliases,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )

    def _resolve_aliases(
        self, genome_digest: str, alias_whitelist: list[str] | None = None
    ) -> list[str]:
        """
        Get the aliases for a genome, optionally restricted to a whitelist.

        Args:
            genome_digest: The digest of the genome.
            alias_whitelist: The list of aliases to include.

        Returns:
            list[str]: The matching alias names.
        """
        return [
            alias
            for alias in self.alias.get_for_genome(genome_digest)
            if (alias_whitelist is None or alias in alias_whitelist)
        ]

    def get_build_paths(
        self,
        genome_digest: str,
        asset_group_name: str | None = None,
        asset_name: str | None = None,
        alias_whitelist: list[str] | None = None,
    ) -> dict[str, Path]:
        """
        Get path to the build directory for the selected genome-group-asset.

        Mirrors :meth:`get_symlink_paths`, but rooted at the ``builds/`` tree.
        Passed as a callback into the removal paths so build bookkeeping is
        torn down with the asset it describes.

        Args:
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            alias_whitelist: The list of aliases to include.

        Returns:
            dict[str, Path]: Mapping of alias names to their build directory paths.
        """
        aliases = self._resolve_aliases(genome_digest, alias_whitelist)
        return get_build_paths_util(
            genome_folder=self.genome_folder,
            aliases=aliases,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )

    def find_build_dir(
        self,
        genome_digest: str,
        asset_group_name: str,
        asset_name: str,
    ) -> Path | None:
        """
        Find the build directory for an asset, whichever alias it was built under.

        The ``builds/`` tree is keyed by the alias used at BUILD time, which is
        not necessarily the alias a later command names. A genome with aliases
        ``hg38`` and ``GRCh38`` built as ``hg38`` has bookkeeping only under
        ``builds/hg38/``, so deriving the path from the alias the user happens
        to type would silently miss it. Search every alias for the genome.

        Args:
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.

        Returns:
            Path | None: The build directory, or None if the asset has no build
            bookkeeping (it was pulled rather than built, or was built before
            this tree existed).
        """
        for path in self.get_build_paths(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        ).values():
            if path.is_dir():
                return path
        return None

    def _symlink_alias(
        self,
        alias_name: str,
        asset_group_name: str | None = None,
        asset_name: str | None = None,
        link_fun: Callable[[str, str], None] = lambda t, s: os.symlink(t, s),  # noqa: E731 - lambda documents (target, symlink) arg order
        alias_whitelist: list[str] | None = None,
    ) -> list[Path]:
        """
        Go through the files in the asset directory and recreate the asset
        directory tree, but instead of copying files, create symbolic links.

        Args:
            alias_name: Alias name for the genome.
            asset_group_name: Asset group name.
            asset_name: Asset name.
            link_fun: Function to use to link files, e.g os.symlink or os.link.
            alias_whitelist: List of aliases to include.

        Returns:
            list[Path]: List of created alias directory paths.
        """
        genome_digest = self.alias.resolve(alias_name)
        if not asset_group_name:
            # Genome-level: render the whole name-addressed alias tree from the
            # assetname rows. The data directory is digest-named, so walking it
            # would produce digest-named alias directories instead of names.
            self.asset.render_alias_tree(genome_digest, alias_whitelist=alias_whitelist)
            return []

        # create symlinks for one asset
        asset_name = asset_name or self.asset.get_default(
            asset_group_name, genome_digest=genome_digest
        )
        src_path = self.asset.get_asset_dir(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )

        target_paths_mapping = self.get_symlink_paths(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            alias_whitelist=alias_whitelist,
        )

        return create_alias_symlinks(
            src_path=src_path,
            target_paths_mapping=target_paths_mapping,
            genome_digest=genome_digest,
            link_fun=link_fun,
        )
