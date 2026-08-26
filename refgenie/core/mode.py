"""
Mode strategy classes for Refgenie.

Encapsulates all behavior that differs between local client and server modes.
The store surface is federated: each mode builds a
:class:`~refgenie.core.store_router.RefgetStoreRouter` holding one or more
:class:`~refgenie.db.tables.Store` backends. Local mode wraps its single on-disk
store as one synthetic ``local`` row; server mode wraps every enabled ``Store``
row from the registry table.
"""

import tempfile
from abc import ABC, abstractmethod
from pathlib import Path
from collections.abc import Callable

from sqlalchemy.engine import Engine

from refgenie.core.store_router import RefgetStoreRouter
from refgenie.db.tables import Store, StoreType
from refgenie.managers.alias import AliasManager, FederatedAliasManager, StoreAliasManager
from refgenie.managers.store import StoreManager


class RefgenieMode(ABC):
    """Encapsulates all behavior that differs between local client and server."""

    @abstractmethod
    def create_store_router(self, store_manager: StoreManager) -> RefgetStoreRouter:
        """Build the federated store router for this mode."""
        ...

    @abstractmethod
    def create_alias_manager(self):
        """Create and return the alias manager for this mode."""
        ...

    @property
    @abstractmethod
    def sequences_enabled(self) -> bool:
        """Whether sequence retrieval (getseq) is available."""
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
        store_path = genome_folder / ".refget_store"
        local = Store(
            name="local",
            url=str(store_path),
            type=StoreType.on_disk,
            priority=0,
        )
        # on_disk stores never touch the cache dir; a genome-folder-local path
        # keeps any stray cache out of the system temp area.
        return RefgetStoreRouter([local], genome_folder / ".refget_store_cache")

    def create_alias_manager(self):
        return FederatedAliasManager(
            StoreAliasManager(genome_folder_getter=self._genome_folder_getter),
            AliasManager(self._database_engine),
        )

    @property
    def sequences_enabled(self) -> bool:
        return True


class ServerMode(RefgenieMode):
    """Server mode. Federated remote RefgetStores. Aliases live in SQL.

    There is no local store here, so the plain SQL manager is the whole alias
    space; a second, permanently empty backend would only add a failure mode.
    """

    def __init__(self, database_engine: Engine):
        self._database_engine = database_engine
        self._cache_dir: Path | None = None

    def create_store_router(self, store_manager: StoreManager) -> RefgetStoreRouter:
        if self._cache_dir is None:
            self._cache_dir = Path(tempfile.mkdtemp(prefix="refgenie_store_cache_"))
        return RefgetStoreRouter(store_manager.enabled_stores(), self._cache_dir)

    def create_alias_manager(self):
        return AliasManager(self._database_engine)

    @property
    def sequences_enabled(self) -> bool:
        return False
