"""Federated store router.

Opens one refget store per enabled :class:`~refgenie.db.tables.Store` row, in
priority order, and dispatches reads to the store that owns a genome. In server
mode (``readonly=True``) each store is loaded in full at startup and frozen, so
no request touches the network; store updates need a redeploy. In local mode
(``readonly=False``) it holds the editable on-disk store.
"""

import time
from pathlib import Path

from refget.store import ReadonlyRefgetStore, RefgetStore

from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError
from refgenie.logger import logger
from refgenie.models import GenomeDigest

#: A store the router holds: editable (local mode) or frozen (server mode).
AnyRefgetStore = RefgetStore | ReadonlyRefgetStore


def make_store(store: Store, cache_dir: Path, *, readonly: bool) -> AnyRefgetStore:
    """Open the refget store for one :class:`Store` row.

    With ``readonly=True`` the store is loaded in full (collections and, for a
    remote store, FHR sidecars) and returned frozen. A failure raises, so a
    server never starts with half a catalog.

    With ``readonly=False``, an ``on_disk`` store gets its collection metadata
    loaded and a ``remote`` store only its collection index; ``store sync``
    pulls the rest per genome.
    """
    t0 = time.monotonic()
    if store.type == StoreType.on_disk:
        rs = RefgetStore.on_disk(str(store.url))
        rs.set_quiet(True)
        rs.load_all_collections()
        if not readonly:
            return rs
        t1 = time.monotonic()
        frozen = rs.into_readonly()
        logger.info(
            f"Store '{store.name}' (on_disk): {frozen.stats().get('n_collections')} collections "
            f"loaded in {t1 - t0:.1f}s"
        )
        return frozen

    cache_dir.mkdir(parents=True, exist_ok=True)
    rs = RefgetStore.open_remote(str(cache_dir), store.url)
    rs.set_quiet(True)
    if not readonly:
        return rs
    t1 = time.monotonic()
    rs.load_all_collections()
    t2 = time.monotonic()
    fhr = rs.pull_fhr()
    t3 = time.monotonic()
    frozen = rs.into_readonly()
    pulled = fhr.get("pulled") if isinstance(fhr, dict) else None
    not_found = fhr.get("not_found") if isinstance(fhr, dict) else None
    logger.info(
        f"Store '{store.name}' (remote {store.url}): "
        f"{frozen.stats().get('n_collections')} collections; "
        f"FHR pulled={pulled} not_found={not_found}; "
        f"open {t1 - t0:.1f}s, collections {t2 - t1:.1f}s, FHR {t3 - t2:.1f}s"
    )
    return frozen


class RefgetStoreRouter:
    """Holds one opened store per enabled :class:`Store` row, priority-ordered."""

    def __init__(self, stores: list[Store], cache_dir: Path, *, readonly: bool):
        self._cache_dir = Path(cache_dir)
        self._stores: dict[str, AnyRefgetStore] = {}
        self._rows: dict[str, Store] = {}
        #: Store names in priority order (lowest ``priority`` integer first).
        self._order: list[str] = []
        for store in sorted(stores, key=lambda s: s.priority):
            logger.debug(f"Opening store '{store.name}' ({store.type.value}: {store.url})")
            self._stores[store.name] = make_store(
                store, self._cache_dir / store.name, readonly=readonly
            )
            self._rows[store.name] = store
            self._order.append(store.name)

    # --- Introspection ------------------------------------------------------

    @property
    def names(self) -> list[str]:
        """Store names in priority order."""
        return list(self._order)

    def get_store(self, name: str) -> AnyRefgetStore:
        """Return the opened store registered under ``name``."""
        if name not in self._stores:
            raise KeyError(f"No store named '{name}' in this router")
        return self._stores[name]

    @property
    def default_store(self) -> AnyRefgetStore:
        """The highest-priority store (the writable one in local mode).

        For a specific genome, use :meth:`store_for_genome`.
        """
        if not self._order:
            raise MissingGenomeError("no stores are configured")
        return self._stores[self._order[0]]

    def iter_stores(self):
        """Yield ``(name, store)`` pairs in priority order."""
        for name in self._order:
            yield name, self._stores[name]

    def store_configs(self) -> list[Store]:
        """The :class:`Store` rows this router opened, in priority order."""
        return [self._rows[name] for name in self._order]

    # --- Dispatch -----------------------------------------------------------

    def store_for_genome(
        self, genome_digest: GenomeDigest, store_name: str | None = None
    ) -> AnyRefgetStore:
        """Return the store that owns ``genome_digest``.

        If ``store_name`` names a known store, it wins outright. Otherwise take
        the first store, in priority order, that has the collection loaded.
        Raises :class:`MissingGenomeError` if no store holds it.
        """
        if store_name and store_name in self._stores:
            return self._stores[store_name]
        for name in self._order:
            if self._stores[name].is_collection_loaded(genome_digest):
                return self._stores[name]
        raise MissingGenomeError(genome_digest)

    # --- Genome-routed pass-throughs ---------------------------------------

    def get_sequence_by_name(
        self, genome_digest: GenomeDigest, name: str, store_name: str | None = None
    ):
        """Fetch a sequence record from the store that owns ``genome_digest``."""
        return self.store_for_genome(genome_digest, store_name).get_sequence_by_name(
            genome_digest, name
        )

    def get_collection_level2(self, genome_digest: GenomeDigest, store_name: str | None = None):
        """Fetch level-2 collection data, or ``None`` if no store owns it."""
        try:
            store = self.store_for_genome(genome_digest, store_name)
        except MissingGenomeError:
            return None
        return store.get_collection_level2(genome_digest)

    def get_fhr_metadata(self, genome_digest: GenomeDigest, store_name: str | None = None):
        """Fetch FHR metadata, or ``None`` if no store owns it or none is set."""
        try:
            store = self.store_for_genome(genome_digest, store_name)
        except MissingGenomeError:
            return None
        return store.get_fhr_metadata(genome_digest)
