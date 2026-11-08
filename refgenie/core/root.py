"""
Refgenie is a SQLModel backed reference genome manager.
"""

import threading
from pathlib import Path

from refget.store import RefgetStore
from refgenie.db.engine import Engine as SqlalchemyDatabaseEngine

from refgenie.config import config
from refgenie.config.settings import LogLevel
from refgenie.const import ALIAS_DIR
from refgenie.core.paths import AssetPaths
from refgenie.core.populate import Populateable, populate_registry_paths
from refgenie.db.events import register_events
from refgenie.db.tables import Alias
from refgenie.exceptions import MissingAliasError, MissingGenomeError
from refgenie.logger import logger
from refgenie.managers import (
    StageManager,
    AssetClassManager,
    AssetManager,
    ConfigurationManager,
    GenomeManager,
    RecipeManager,
    RemoteManager,
)
from refgenie.managers.alias import AliasBackend
from refgenie.db.engine import create_database_engine
from refgenie.managers.database import DatabaseManager
from refgenie.managers.sequence import SequenceManager
from refgenie.managers.transfer import TransferManager
from refgenie.managers.build import BuildManager
from refgenie.managers.transfer.puller import AssetPuller
from refgenie.managers.plugin_settings import PluginSettingsManager
from refgenie.managers.species import SpeciesManager
from refgenie.managers.store import StoreManager
from refgenie.managers.sources import RefgenieserverClient, ServerManager, SourceManager
from refgenie.core.mode import LocalMode, RefgenieMode, ServerMode
from refgenie.core.store_router import RefgetStoreRouter
from refgenie.plugins.events import EventSink, update_scope
from refgenie.plugins.host import PluginHost
from refgenie.models import (
    AssetRegistryPathComponents,
    GenomeAlias,
    GenomeDigest,
)


class Refgenie:
    """
    A reference genome manager.

    Callers use the managers directly (``rgc.alias``, ``rgc.asset``,
    ``rgc.transfer``, ...). Methods here only coordinate several managers; see
    "AssetManager is composed, not mixed" in ``docs/design-notes.md``.
    """

    def __init__(
        self,
        database_config_path: str | Path | None = None,
        database_engine: SqlalchemyDatabaseEngine | None = None,
        server_clients_mapping: dict[str, RefgenieserverClient] | None = None,
        suppress_migrations: bool = False,
        server_mode: bool = False,
        plugins: bool | None = None,
    ):
        """
        Initialize the reference genome manager.

        Args:
            database_config_path: The database configuration file; defaults to
                the configured path.
            database_engine: The database engine; defaults to one built from the config.
            server_clients_mapping: Server URL -> client, handed to ``rgc.servers``.
                Clients are otherwise created on first use.
            suppress_migrations: Do not apply migrations even if needed.
            server_mode: Run in server mode (federated stores from the ``store``
                table, SQL aliases, no sequence ingestion) instead of local mode.
            plugins: Whether to run installed plugins. None means on in local
                mode and off in server mode unless ``REFGENIE_SERVER_PLUGINS``
                is set. ``REFGENIE_DISABLE_PLUGINS`` always turns them off.
        """
        register_events()
        logger.debug(f"{config=}")
        self._database_engine = database_engine or create_database_engine(
            database_config_path=database_config_path,
            echo=config.log_level == LogLevel.DEBUG,
        )
        self._server_clients_mapping = server_clients_mapping
        self.recipe = RecipeManager(self.database_engine)
        self.asset_class = AssetClassManager(self.database_engine)
        self.configuration = ConfigurationManager(self.database_engine)
        #: Schema creation, backend init, migrations, and purge.
        self.database = DatabaseManager(self.database_engine, self.configuration)
        #: Push remotes and which staged assets go to each.
        self.remote = RemoteManager(self.database_engine)
        self.stage = StageManager(self.database_engine, self.remote)
        #: Data channels, which publish recipes and asset classes, and syncing them.
        self.sources = SourceManager(self.database_engine, self.asset_class, self.recipe)
        # PluginHost does no discovery, so this costs nothing until a hook fires.
        if plugins is None:
            plugins = config.server_plugins if server_mode else True
        self._plugins_enabled = plugins
        self._events = EventSink(PluginHost(self) if plugins else None)
        self._plugin_settings_manager: PluginSettingsManager | None = None
        if server_mode:
            self._mode: RefgenieMode = ServerMode(
                genome_folder_getter=lambda: self.genome_folder,
                database_engine=self._database_engine,
            )
        else:
            self._mode = LocalMode(
                genome_folder_getter=lambda: self.genome_folder,
                database_engine=self._database_engine,
            )

        #: The federation registry of RefgetStores.
        self.store = StoreManager(self.database_engine)
        #: Species names from the stores' species.json, for search.
        self.species = SpeciesManager(self.database_engine)
        # Lazy: server mode picks its alias manager by whether the genome folder
        # holds a store, and the genome folder is only known after migrations.
        self._alias_manager: AliasBackend | None = None
        self._store_router: RefgetStoreRouter | None = None  # Lazy, built by mode
        # Server startup warms the router on a background thread while requests
        # may already be asking for it; build it once.
        self._store_router_lock = threading.Lock()

        self.genome = GenomeManager(
            self.database_engine,
            refget_store_getter=lambda: self.refget_store,
            alias_manager_getter=lambda: self.alias,
            servers_getter=lambda: self.servers,
            events=self._events,
            store_router_getter=lambda: self.store_router,
            species_manager=self.species,
        )
        #: Sequence retrieval through the store router.
        self.sequence = SequenceManager(
            genome_manager=self.genome,
            store_router_getter=lambda: self.store_router,
            sequences_enabled=self._mode.sequences_enabled,
        )
        # Lazy: each needs the genome folder, which is read from the database
        # and so is only available after migrations.
        self._asset_manager: AssetManager | None = None
        self._server_manager: ServerManager | None = None
        self._build_manager: BuildManager | None = None
        self._transfer_manager: TransferManager | None = None

        if not suppress_migrations:
            self.database.migrate()

    def __str__(self):
        clients = [] if self._server_manager is None else list(self._server_manager.clients)
        return f"Refgenie(database_engine={self.database_engine}, server_clients={clients})"

    def __repr__(self):
        return str(self)

    @property
    def store_router(self) -> RefgetStoreRouter:
        """The federated store router, built by the mode on first access."""
        if self._store_router is None:
            with self._store_router_lock:
                if self._store_router is None:
                    self._store_router = self._mode.create_store_router(self.store)
        return self._store_router

    def reload_store_router(self) -> RefgetStoreRouter:
        """Rebuild the router from the current ``store`` registry.

        Call after ``store add``/``sync``/``remove``. The alias manager caches
        against the store it reads, so it is invalidated too.

        The new router is built before it replaces the old one, so a server's
        concurrent requests keep reading the old router meanwhile, and a store
        that fails to load leaves the old router in place.
        """
        router = self._mode.create_store_router(self.store)
        with self._store_router_lock:
            self._store_router = router
        self.alias.invalidate()
        return router

    @property
    def refget_store(self) -> RefgetStore:
        """The default (highest-priority, writable) store.

        Writes go here; genome-specific reads should go through :attr:`store_router`.
        """
        return self.store_router.default_store

    @property
    def alias(self) -> AliasBackend:
        """The alias manager the mode selects."""
        if self._alias_manager is None:
            self._alias_manager = self._mode.create_alias_manager(
                lambda: self.refget_store, self._events
            )
        return self._alias_manager

    @property
    def plugins(self) -> PluginSettingsManager:
        """Plugin settings, and what plugins are installed. See ``refgenie/plugins/``."""
        if self._plugin_settings_manager is None:
            self._plugin_settings_manager = PluginSettingsManager(
                self.database_engine, enabled=self._plugins_enabled
            )
        return self._plugin_settings_manager

    def batch_updates(self):
        """Group every change made inside this block into one ``post_update``.

        Each marked operation already fires its own ``post_update``; use this
        when one user action runs several of them (the CLI wraps every command).
        """
        return self._events.scope()

    @property
    def asset(self) -> AssetManager:
        """The asset manager: the public interface for all asset operations."""
        if self._asset_manager is None:
            self._asset_manager = AssetManager(
                database_engine=self.database_engine,
                genome_folder=self.genome_folder,
                alias_folder=self.alias_folder,
                alias_manager=self.alias,
                genome_manager=self.genome,
                asset_class_manager=self.asset_class,
                events=self._events,
            )
        return self._asset_manager

    def paths(self) -> AssetPaths:
        """A fresh lazy view of every local asset path. See ``core/paths.py``.

        A method, not a cached property: each caller owns the view's memo, so a
        long-running process never serves paths from before a pull or remove.
        """
        return AssetPaths(self)

    @property
    def servers(self) -> ServerManager:
        """The refgenieservers this node pulls from: subscriptions, clients, catalog, remote seek."""
        if self._server_manager is None:
            self._server_manager = ServerManager(
                database_engine=self.database_engine,
                alias_manager=self.alias,
                seek_keys=self.asset.seek_key,
                clients=self._server_clients_mapping,
            )
        return self._server_manager

    @property
    def build(self) -> BuildManager:
        """Building assets from recipes: build, preflight, and Snakemake targets."""
        if self._build_manager is None:
            self._build_manager = BuildManager(
                database_engine=self.database_engine,
                recipe_manager=self.recipe,
                asset=self.asset,
                alias_manager=self.alias,
                genome_manager=self.genome,
                sequence_manager=self.sequence,
                stage_manager=self.stage,
                stage_folder_getter=lambda: self.genome_stage_folder,
                pull_parent=lambda alias: self.transfer.pull(
                    asset_group_name="fasta", genome=alias
                ),
                events=self._events,
            )
        return self._build_manager

    @property
    def transfer(self) -> TransferManager:
        """Pulls from subscribed servers: one asset, many genomes, or a mirror."""
        if self._transfer_manager is None:
            self._transfer_manager = TransferManager(
                servers=self.servers,
                alias_manager=self.alias,
                genome_manager=self.genome,
                puller=AssetPuller(
                    database_engine=self.database_engine,
                    asset=self.asset,
                    servers=self.servers,
                    alias_manager=self.alias,
                    genome_manager=self.genome,
                    events=self._events,
                    asset_class=self.asset_class,
                ),
                events=self._events,
            )
        return self._transfer_manager

    @property
    def database_engine(self) -> SqlalchemyDatabaseEngine:
        """The database engine."""
        return self._database_engine

    @property
    def genome_folder(self) -> Path:
        """The genome folder, from the latest configuration."""
        configuration = self.configuration.get_latest()
        return Path(configuration.genome_folder)

    @property
    def genome_stage_folder(self) -> Path | None:
        """The genome stage folder, or None if not set."""
        configuration = self.configuration.get_latest()
        return None if (a := configuration.genome_stage_folder) is None else Path(a)

    @property
    def data_folder(self) -> Path:
        """The data folder inside the genome folder."""
        return self.genome_folder / "data"

    @property
    def alias_folder(self) -> Path:
        """The alias folder inside the genome folder."""
        return self.genome_folder / ALIAS_DIR

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

    def populate(self, input: Populateable):
        """Replace ``refgenie://`` registry paths in the input with local paths.

        Args:
            input: A string, or a list or dict of them, nested to any depth.

        Returns:
            The input, with each registry path replaced.
        """
        return populate_registry_paths(self, input)

    def populater(self, input: Populateable, server_urls: list[str] | None = None):
        """Replace ``refgenie://`` registry paths in the input with remote URLs.

        Args:
            input: A string, or a list or dict of them, nested to any depth.
            server_urls: Servers to ask. Defaults to the subscriptions.

        Returns:
            The input, with each registry path replaced.
        """
        return populate_registry_paths(self, input, remote=True, server_urls=server_urls)

    @update_scope
    def set_genome_alias(
        self,
        alias_name: GenomeAlias,
        genome_digest: GenomeDigest | None = None,
        genome_description: str | None = None,
        server_urls: list[str] | None = None,
    ) -> Alias:
        """
        Set a genome alias, asking the servers for the digest if none is given.

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
        # Make sure the genome exists locally. The alias is added below, so
        # pass no alias_names to the genome layer.
        if genome_digest is None:
            genome_digest = self.servers.resolve_alias(alias_name, server_urls)
            if genome_digest is None:
                raise MissingAliasError(alias_name)
            logger.info(f"Determined digest for {alias_name}: {genome_digest}")
            source = self.servers.genome_source(server_urls)
            if source is None:
                raise ConnectionError(f"No server can serve as a genome source for '{alias_name}'.")

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

        alias = self.alias.add(alias_name, genome_digest)
        self.asset.tree.render_genome(genome_digest, aliases=[alias_name])
        return alias
