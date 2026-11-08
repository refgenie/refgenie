"""StoreManager - CRUD for the ``store`` federation registry.

The ``store`` table lists the refget stores a refgenie service federates over.
The root reads its enabled rows to build the
:class:`~refgenie.core.store_router.RefgetStoreRouter`.
"""

from collections.abc import Iterable

from rich.table import Table

from refgenie.db.engine import Engine
from refgenie.db.queries import alias as alias_queries
from refgenie.db.queries import genome as genome_queries
from refgenie.db.queries import store as queries
from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import StoreExistsError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
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
            if queries.find(session, name) is not None:
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
            return queries.get(session, name)

    def remove(self, name: str) -> None:
        """Remove a store row. Genomes keep their (now dangling) ``store_name``."""
        with self._database_session as session:
            store = queries.get(session, name)
            session.delete(store)
            session.commit()
            logger.info(f"Removed store '{name}'")

    def list_all(self) -> Iterable[Store]:
        """All store rows, priority order."""
        with self._database_session as session:
            return queries.list_all(session)

    def enabled_stores(self) -> list[Store]:
        """Enabled store rows in priority order (lowest integer first)."""
        with self._database_session as session:
            return queries.list_all(session, enabled_only=True)

    def genome_counts(self) -> dict[str, int]:
        """Map store name -> number of genomes it owns."""
        with self._database_session as session:
            return genome_queries.count_by_store(session)

    def apply_declaration(
        self, declared: list[dict], disable_missing: bool = False
    ) -> dict[str, list[str]]:
        """
        Idempotently upsert each declared store into the registry.

        Insert if new; update url/type/priority/enabled/description when they
        drift. Rows are never deleted (genomes keep a foreign key to them).
        With ``disable_missing``, every existing store not in the declaration
        is set ``enabled=False``. An entry with no name or url, or with the
        reserved name, is logged and skipped.

        Args:
            declared: Store objects, each ``{"name", "url", "type"?,
                "priority"?, "description"?, "enabled"?}``.
            disable_missing: Disable stores the declaration does not list.

        Returns:
            Store names by change: ``{"added", "updated", "disabled"}``.
        """
        changes: dict[str, list[str]] = {"added": [], "updated": [], "disabled": []}
        declared_names: set[str] = set()
        with self._database_session as session:
            for entry in declared:
                name = entry.get("name")
                url = entry.get("url")
                if not name or not url:
                    logger.error(f"Skipping store declaration with no name/url: {entry!r}")
                    continue
                try:
                    check_store_name(name)
                except ValueError as e:
                    logger.error(f"Skipping store declaration: {e}")
                    continue
                declared_names.add(name)
                values = {
                    "url": url,
                    "type": StoreType(entry.get("type", "remote")),
                    "priority": int(entry.get("priority", 100)),
                    "enabled": bool(entry.get("enabled", True)),
                    "description": entry.get("description"),
                }

                existing = queries.find(session, name)
                if existing is None:
                    session.add(Store(name=name, **values))
                    changes["added"].append(name)
                    logger.info(
                        f"Registered store '{name}' "
                        f"({values['type'].value}, priority={values['priority']})"
                    )
                elif any(getattr(existing, k) != v for k, v in values.items()):
                    for k, v in values.items():
                        setattr(existing, k, v)
                    changes["updated"].append(name)
                    logger.info(f"Updated store '{name}'")
            if disable_missing:
                for store in queries.list_all(session, enabled_only=True):
                    if store.name not in declared_names:
                        store.enabled = False
                        changes["disabled"].append(store.name)
                        logger.info(f"Disabled store '{store.name}' (not declared)")
            session.commit()
        return changes

    def conflicts(self) -> list[tuple[str, str, str]]:
        """
        Alias collisions: the losers stored under qualified ``store::alias`` names.

        Returns:
            ``(alias, store, genome digest)`` for each loser, by qualified name.
        """
        with self._database_session as session:
            aliases = alias_queries.list_all(session)
        rows = []
        for alias in aliases:
            if "::" in alias.name:
                store_name, _, bare = alias.name.partition("::")
                rows.append((bare, store_name, alias.genome_digest))
        return rows

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
