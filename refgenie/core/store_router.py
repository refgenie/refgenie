"""Federated store router.

Refgenie serves genomes from several physically separate refget stores as one
service (Design B). This module opens one :class:`~refget.store.RefgetStore` per
enabled :class:`~refgenie.db.tables.Store` row, in priority order, and dispatches
reads to the store that owns a given genome.

Content is digest-addressed, so identical genomes across stores dedup for free.
The only true conflict is two stores mapping the same alias to different
collection digests, resolved by store priority (lower integer wins); that policy
lives in the alias layer, not here.
"""

from pathlib import Path

from refget.store import RefgetStore

from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError
from refgenie.logger import logger


def make_store(store: Store, cache_dir: Path) -> RefgetStore:
    """Open the :class:`RefgetStore` for one :class:`Store` row.

    ``on_disk`` stores are opened directly and only their collection metadata is
    loaded (never sequence bytes -- listing and alias resolution never need
    them, and a metadata-only store has no bytes on disk). ``remote`` stores are
    opened over http(s)/S3 and additionally pull their alias and FHR sidecars,
    exactly as the old single-store server path did.
    """
    if store.type == StoreType.on_disk:
        rs = RefgetStore.on_disk(str(store.url))
        rs.set_quiet(True)
        rs.load_all_collections()
        return rs

    cache_dir.mkdir(parents=True, exist_ok=True)
    rs = RefgetStore.open_remote(str(cache_dir), store.url)
    rs.set_quiet(True)
    rs.load_all_collections()
    try:
        rs.pull_aliases()
    except Exception as e:  # noqa: BLE001 - a store with no aliases is tolerated
        logger.debug(f"pull_aliases failed for store '{store.name}' ({store.url}): {e}")
    try:
        rs.pull_fhr()
    except Exception as e:  # noqa: BLE001 - a store with no FHR is tolerated
        logger.debug(f"pull_fhr failed for store '{store.name}' ({store.url}): {e}")
    return rs


class RefgetStoreRouter:
    """Holds one opened :class:`RefgetStore` per enabled store, priority-ordered.

    A single-store deployment (local mode) is just a router with one store; the
    federated server mode wraps every enabled :class:`Store` row.
    """

    def __init__(self, stores: list[Store], cache_dir: Path):
        self._cache_dir = Path(cache_dir)
        self._stores: dict[str, RefgetStore] = {}
        self._rows: dict[str, Store] = {}
        #: Store names in priority order (lowest ``priority`` integer first).
        self._order: list[str] = []
        for store in sorted(stores, key=lambda s: s.priority):
            logger.debug(f"Opening store '{store.name}' ({store.type.value}: {store.url})")
            self._stores[store.name] = make_store(store, self._cache_dir / store.name)
            self._rows[store.name] = store
            self._order.append(store.name)

    # --- Introspection ------------------------------------------------------

    @property
    def names(self) -> list[str]:
        """Store names in priority order."""
        return list(self._order)

    def get_store(self, name: str) -> RefgetStore:
        """Return the opened store registered under ``name``."""
        if name not in self._stores:
            raise KeyError(f"No store named '{name}' in this router")
        return self._stores[name]

    @property
    def default_store(self) -> RefgetStore:
        """The highest-priority store.

        This is the writable store in local mode (the single on-disk store) and
        the default mount / write-path target elsewhere. Callers that operate on
        one specific genome should prefer :meth:`store_for_genome`.
        """
        if not self._order:
            raise MissingGenomeError("no stores are configured")
        return self._stores[self._order[0]]

    def iter_stores(self):
        """Yield ``(name, RefgetStore)`` pairs in priority order."""
        for name in self._order:
            yield name, self._stores[name]

    def store_configs(self) -> list[Store]:
        """The :class:`Store` rows this router opened, in priority order."""
        return [self._rows[name] for name in self._order]

    # --- Dispatch -----------------------------------------------------------

    def store_for_genome(self, genome_digest: str, store_name: str | None = None) -> RefgetStore:
        """Return the store that owns ``genome_digest``.

        If ``store_name`` names a known store, it wins outright. Otherwise probe
        each store in priority order and return the first that has the collection
        loaded. Raises :class:`MissingGenomeError` if no store holds it.
        """
        if store_name and store_name in self._stores:
            return self._stores[store_name]
        for name in self._order:
            store = self._stores[name]
            try:
                if store.is_collection_loaded(genome_digest):
                    return store
            except Exception as e:  # noqa: BLE001 - a probe failure is not fatal
                logger.debug(f"is_collection_loaded probe failed on store '{name}': {e}")
        raise MissingGenomeError(genome_digest)

    # --- Genome-routed pass-throughs ---------------------------------------

    def get_sequence_by_name(self, genome_digest: str, name: str, store_name: str | None = None):
        """Fetch a sequence record from the store that owns ``genome_digest``."""
        return self.store_for_genome(genome_digest, store_name).get_sequence_by_name(
            genome_digest, name
        )

    def get_collection_level2(self, genome_digest: str, store_name: str | None = None):
        """Fetch level-2 collection data, or ``None`` if no store owns it."""
        try:
            store = self.store_for_genome(genome_digest, store_name)
        except MissingGenomeError:
            return None
        return store.get_collection_level2(genome_digest)

    def get_fhr_metadata(self, genome_digest: str, store_name: str | None = None):
        """Fetch FHR metadata, or ``None`` if no store owns it or none is set."""
        try:
            store = self.store_for_genome(genome_digest, store_name)
        except MissingGenomeError:
            return None
        return store.get_fhr_metadata(genome_digest)
