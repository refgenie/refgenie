"""StoreManager - CRUD for the ``store`` federation registry.

The ``store`` table is the single source of truth for which refget stores a
refgenie service federates over. ``refgenie store add/sync/list/remove`` drive
this manager, and the root reads its enabled rows to build the
:class:`~refgenie.core.store_router.RefgetStoreRouter`.
"""

from collections.abc import Iterable

from rich.table import Table
from sqlalchemy import func
from sqlalchemy.engine import Engine
from sqlmodel import select

from refgenie.db.tables import Genome, Store, StoreType
from refgenie.exceptions import MissingStoreError, StoreExistsError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.utils.tables import build_table

#: Router name of the store a node builds into (``<genome_folder>/.refget_store``).
#: Reserved: server mode federates that store under this name, so a registered
#: store may not take it.
LOCAL_STORE_NAME = "local"


def check_store_name(name: str) -> None:
    """Raise ``ValueError`` if ``name`` is reserved for the node's own store."""
    if name == LOCAL_STORE_NAME:
        raise ValueError(
            f"Store name '{LOCAL_STORE_NAME}' is reserved for the genome folder's own store."
        )


class StoreManager(ResourceManager):
    """Manager for the ``store`` federation registry."""

    def __init__(self, database_engine: Engine):
        super().__init__(database_engine)

    def add(
        self,
        name: str,
        url: str,
        store_type: StoreType = StoreType.remote,
        priority: int = 100,
        enabled: bool = True,
        description: str | None = None,
    ) -> Store:
        """Register a store.

        Raises :class:`StoreExistsError` on a name clash, and ``ValueError`` for
        the reserved name :data:`LOCAL_STORE_NAME`.
        """
        check_store_name(name)
        with self._database_session as session:
            existing = session.exec(select(Store).where(Store.name == name)).one_or_none()
            if existing is not None:
                raise StoreExistsError(name)
            store = Store(
                name=name,
                url=url,
                type=store_type,
                priority=priority,
                enabled=enabled,
                description=description,
            )
            session.add(store)
            session.commit()
            session.refresh(store)
            logger.info(f"Added store '{name}' ({store_type.value}, priority={priority})")
            return store

    def get(self, name: str) -> Store:
        """Return the store row named ``name``."""
        with self._database_session as session:
            return one_or_raise(
                session,
                select(Store).where(Store.name == name),
                MissingStoreError(name),
                unique=True,
            )

    def remove(self, name: str) -> None:
        """Remove a store row. Genomes keep their (now dangling) ``store_name``."""
        with self._database_session as session:
            store = one_or_raise(
                session,
                select(Store).where(Store.name == name),
                MissingStoreError(name),
                unique=True,
            )
            session.delete(store)
            session.commit()
            logger.info(f"Removed store '{name}'")

    def list_all(self) -> Iterable[Store]:
        """All store rows, priority order."""
        with self._database_session as session:
            return session.exec(select(Store).order_by(Store.priority)).all()

    def enabled_stores(self) -> list[Store]:
        """Enabled store rows in priority order (lowest integer first)."""
        with self._database_session as session:
            return list(
                session.exec(select(Store).where(Store.enabled).order_by(Store.priority)).all()
            )

    def genome_counts(self) -> dict[str, int]:
        """Map store name -> number of genomes it owns."""
        with self._database_session as session:
            rows = session.exec(
                select(Genome.store_name, func.count(Genome.digest)).group_by(Genome.store_name)
            ).all()
        return {name: count for name, count in rows if name is not None}

    def table(self) -> Table:
        """A Rich table of all stores with their genome counts."""
        counts = self.genome_counts()
        rows = []
        for store in self.list_all():
            rows.append(
                (
                    store.name,
                    store.url,
                    store.type.value,
                    str(store.priority),
                    "yes" if store.enabled else "no",
                    str(counts.get(store.name, 0)),
                )
            )
        return build_table(
            "Stores",
            ["Name", "URL", "Type", "Priority", "Enabled", "Genomes"],
            rows,
        )
