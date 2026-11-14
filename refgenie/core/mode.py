"""
Mode strategy classes: the behavior that differs between local and server mode.

Local mode routes to its single on-disk store, as a synthetic ``local`` row.
Server mode routes to every enabled ``Store`` row, plus that ``local`` row when
the served genome folder holds a store of its own.
"""

import tempfile
from abc import ABC, abstractmethod
from collections.abc import Callable
from pathlib import Path

from refget.store import RefgetStore

from refgenie.core.store_router import RefgetStoreRouter
from refgenie.db.engine import Engine
from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingConfigDataError
from refgenie.managers.alias import (
    AliasBackend,
    AliasManager,
    FederatedAliasManager,
    _StoreAliases,
)
from refgenie.managers.store import LOCAL_STORE_NAME, StoreManager
from refgenie.plugins.events import EventSink


def _local_store_row(genome_folder: Path, priority: int = 0) -> Store:
    """The synthetic registry row for a genome folder's own on-disk store."""
    return Store(
        name=LOCAL_STORE_NAME,
        url=str(genome_folder / ".refget_store"),
        type=StoreType.on_disk,
        priority=priority,
    )


def _federated_aliases(
    refget_store_getter: Callable[[], RefgetStore],
    genome_folder_getter: Callable[[], Path],
    database_engine: Engine,
    events: EventSink | None,
) -> FederatedAliasManager:
    """The local store's aliases, with the SQL alias table behind them."""
    return FederatedAliasManager(
        _StoreAliases(refget_store_getter, genome_folder_getter),
        AliasManager(database_engine),
        events=events,
    )


class RefgenieMode(ABC):
    """Encapsulates all behavior that differs between local client and server."""

    @abstractmethod
    def create_store_router(self, store_manager: StoreManager) -> RefgetStoreRouter:
        """Build the federated store router for this mode."""
        ...

    @abstractmethod
    def create_alias_manager(
        self,
        refget_store_getter: Callable[[], RefgetStore],
        events: EventSink | None = None,
    ) -> AliasBackend:
        """
        Create the alias manager for this mode, fully wired.

        Args:
            refget_store_getter: Returns the default store of the router. A
                getter, because the router is built lazily.
            events: Where alias changes are recorded for plugins. Handed to
                the outermost backend only, so one change is recorded once.
        """
        ...

    @property
    @abstractmethod
    def sequences_enabled(self) -> bool:
        """Whether sequence retrieval (``rgc.sequence.get``) is available."""
        ...


class LocalMode(RefgenieMode):
    """Local client mode. Owns a single editable on-disk RefgetStore.

    Aliases for the genomes this node built live in that store. Aliases for the
    genomes it federates over live in the SQL ``alias`` table, so its alias
    manager reads both.
    """

    def __init__(self, genome_folder_getter: Callable[[], Path], database_engine: Engine):
        self._genome_folder_getter = genome_folder_getter
        self._database_engine = database_engine

    def create_store_router(self, store_manager: StoreManager) -> RefgetStoreRouter:
        genome_folder = self._genome_folder_getter()
        # on_disk stores never touch the cache dir; a genome-folder-local path
        # keeps any stray cache out of the system temp area.
        return RefgetStoreRouter(
            [_local_store_row(genome_folder)],
            genome_folder / ".refget_store_cache",
            readonly=False,
            # Every genome this node reads is in its one store, pulled or built.
            sole_store=True,
        )

    def create_alias_manager(
        self,
        refget_store_getter: Callable[[], RefgetStore],
        events: EventSink | None = None,
    ) -> AliasBackend:
        return _federated_aliases(
            refget_store_getter, self._genome_folder_getter, self._database_engine, events
        )

    @property
    def sequences_enabled(self) -> bool:
        return True


class ServerMode(RefgenieMode):
    """Server mode. Federated read-only stores, SQL aliases, no sequence ingestion.

    If the served genome folder has its own ``.refget_store`` (genomes built on
    this node), it joins as the highest-priority store and its aliases are read
    first, as in local mode. Otherwise aliases come from SQL alone.
    """

    def __init__(self, genome_folder_getter: Callable[[], Path], database_engine: Engine):
        self._genome_folder_getter = genome_folder_getter
        self._database_engine = database_engine
        self._cache_dir: Path | None = None

    def _own_store_folder(self) -> Path | None:
        """The served genome folder, if it holds a store built on this node."""
        try:
            genome_folder = self._genome_folder_getter()
        except MissingConfigDataError:
            return None
        return genome_folder if (genome_folder / ".refget_store").is_dir() else None

    def create_store_router(self, store_manager: StoreManager) -> RefgetStoreRouter:
        if self._cache_dir is None:
            self._cache_dir = Path(tempfile.mkdtemp(prefix="refgenie_store_cache_"))
        stores = store_manager.enabled_stores()
        genome_folder = self._own_store_folder()
        if genome_folder is not None:
            # Rank it first: the alias manager reads the router's default store.
            top = min((s.priority for s in stores), default=1)
            stores = [_local_store_row(genome_folder, priority=top - 1), *stores]
        # Server mode never writes to a store: load each in full, once, frozen.
        return RefgetStoreRouter(stores, self._cache_dir, readonly=True)

    def create_alias_manager(
        self,
        refget_store_getter: Callable[[], RefgetStore],
        events: EventSink | None = None,
    ) -> AliasBackend:
        if self._own_store_folder() is None:
            return AliasManager(self._database_engine, events=events)
        return _federated_aliases(
            refget_store_getter, self._genome_folder_getter, self._database_engine, events
        )

    @property
    def sequences_enabled(self) -> bool:
        return False
