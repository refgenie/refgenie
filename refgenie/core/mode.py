"""
Mode strategy classes for Refgenie.

Encapsulates all behavior that differs between local client and server modes.
The store surface is federated: each mode builds a
:class:`~refgenie.core.store_router.RefgetStoreRouter` holding one or more
:class:`~refgenie.db.tables.Store` backends. Local mode wraps its single on-disk
store as one synthetic ``local`` row; server mode wraps every enabled ``Store``
row from the registry table, plus that same ``local`` row when the served genome
folder holds a store of its own.
"""

import tempfile
from abc import ABC, abstractmethod
from pathlib import Path
from collections.abc import Callable

from refget.store import RefgetStore
from sqlalchemy.engine import Engine

from refgenie.core.store_router import RefgetStoreRouter
from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingConfigDataError
from refgenie.managers.alias import (
    AliasBackend,
    AliasManager,
    FederatedAliasManager,
    StoreAliasManager,
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
    """Local client mode. Owns a single on-disk RefgetStore.

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
            [_local_store_row(genome_folder)], genome_folder / ".refget_store_cache"
        )

    def create_alias_manager(
        self,
        refget_store_getter: Callable[[], RefgetStore],
        events: EventSink | None = None,
    ) -> AliasBackend:
        return FederatedAliasManager(
            StoreAliasManager(refget_store_getter, self._genome_folder_getter),
            AliasManager(self._database_engine),
            events=events,
        )

    @property
    def sequences_enabled(self) -> bool:
        return True


class ServerMode(RefgenieMode):
    """Server mode. Federated RefgetStores, SQL aliases, no sequence ingestion.

    The stores come from the registry table. A served database that was built
    on this node also has its own ``.refget_store``, which holds the collections
    and aliases of every genome built here. That store joins the federation as
    the highest-priority member, and its aliases are read first, the same way
    local mode reads them. A server with no store of its own reads the SQL
    ``alias`` table alone.
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
            # Rank it ahead of every registered store: a genome built on this
            # node owns its name, and the router's default store must be this
            # one, because that is the store the alias manager reads.
            top = min((s.priority for s in stores), default=1)
            stores = [_local_store_row(genome_folder, priority=top - 1), *stores]
        return RefgetStoreRouter(stores, self._cache_dir)

    def create_alias_manager(
        self,
        refget_store_getter: Callable[[], RefgetStore],
        events: EventSink | None = None,
    ) -> AliasBackend:
        if self._own_store_folder() is None:
            return AliasManager(self._database_engine, events=events)
        return FederatedAliasManager(
            StoreAliasManager(refget_store_getter, self._genome_folder_getter),
            AliasManager(self._database_engine),
            events=events,
        )

    @property
    def sequences_enabled(self) -> bool:
        return False
