"""Federated store router.

Refgenie serves genomes from several physically separate refget stores as one
service (Design B). This module opens one :class:`~refget.store.RefgetStore` per
enabled :class:`~refgenie.db.tables.Store` row, in priority order, and dispatches
reads to the store that owns a given genome.

Content is digest-addressed, so identical genomes across stores dedup for free.
The only true conflict is two stores mapping the same alias to different
collection digests, resolved by store priority (lower integer wins); that policy
lives in the alias layer, not here.

Remote stores open lazily: only the collection index (stubs) is fetched up
front. A genome's full collection and its FHR sidecar
are fetched from the remote store the first time a request routes to it
(:meth:`RefgetStoreRouter.store_for_genome`), then served from memory.
"""

import threading
from pathlib import Path

from refget.store import RefgetStore

from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError
from refgenie.logger import logger
from refgenie.models import GenomeDigest


def make_store(store: Store, cache_dir: Path) -> RefgetStore:
    """Open the :class:`RefgetStore` for one :class:`Store` row.

    ``on_disk`` stores are opened directly and only their collection metadata is
    loaded (never sequence bytes -- listing and alias resolution never need
    them, and a metadata-only store has no bytes on disk). ``remote`` stores are
    opened over http(s)/S3 with only their collection index (stubs). Full
    collections and FHR sidecars are fetched one genome at a time on first use
    (see :meth:`RefgetStoreRouter.store_for_genome`): loading them all up front
    costs one S3 request per collection, minutes for a large store.

    Alias sidecars are not pulled either: the server reads aliases from SQL (or
    its own on-disk store), never from a remote store, and pulling them costs
    15-25 s per store. ``store sync``, the one reader, pulls them itself.
    """
    if store.type == StoreType.on_disk:
        rs = RefgetStore.on_disk(str(store.url))
        rs.set_quiet(True)
        rs.load_all_collections()
        return rs

    cache_dir.mkdir(parents=True, exist_ok=True)
    rs = RefgetStore.open_remote(str(cache_dir), store.url)
    rs.set_quiet(True)
    return rs


class RefgetStoreRouter:
    """Holds one opened :class:`RefgetStore` per enabled store, priority-ordered.

    A single-store deployment (local mode) is just a router with one store; the
    federated server mode wraps every enabled :class:`Store` row.

    Remote stores are opened with stubs only. :meth:`store_for_genome` loads the
    requested genome's full collection and FHR sidecar on first use, once per
    ``(store, digest)``, under a per-store lock so concurrent first requests
    (FastAPI runs sync handlers in a threadpool) neither double-load nor race
    inside the Rust store.
    """

    def __init__(self, stores: list[Store], cache_dir: Path):
        self._cache_dir = Path(cache_dir)
        self._stores: dict[str, RefgetStore] = {}
        self._rows: dict[str, Store] = {}
        #: Store names in priority order (lowest ``priority`` integer first).
        self._order: list[str] = []
        #: ``(store_name, digest)`` pairs whose collection and FHR are loaded.
        self._loaded: set[tuple[str, str]] = set()
        self._locks: dict[str, threading.Lock] = {}
        for store in sorted(stores, key=lambda s: s.priority):
            logger.debug(f"Opening store '{store.name}' ({store.type.value}: {store.url})")
            self._stores[store.name] = make_store(store, self._cache_dir / store.name)
            self._rows[store.name] = store
            self._order.append(store.name)
            self._locks[store.name] = threading.Lock()

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

    def _is_remote(self, name: str) -> bool:
        return self._rows[name].type != StoreType.on_disk

    def _holds(self, name: str, genome_digest: GenomeDigest) -> bool:
        """Whether store ``name`` holds ``genome_digest``, without loading it.

        Remote stores hold stubs, for which ``is_collection_loaded`` is False,
        so they are probed through the stub-level collection metadata instead.
        """
        store = self._stores[name]
        if self._is_remote(name):
            return store.get_collection_metadata(genome_digest) is not None
        return store.is_collection_loaded(genome_digest)

    def _ensure_loaded(self, name: str, genome_digest: GenomeDigest) -> bool:
        """Load one collection and its FHR sidecar from a remote store, once.

        Returns False when the store does not hold the digest. A network error
        while loading the collection propagates; a missing FHR sidecar is normal
        and only logged.
        """
        if not self._is_remote(name) or (name, genome_digest) in self._loaded:
            return True
        store = self._stores[name]
        with self._locks[name]:
            if (name, genome_digest) in self._loaded:
                return True
            if store.get_collection_metadata(genome_digest) is None:
                return False
            try:
                store.load_collection(genome_digest)
            except Exception as e:
                raise RuntimeError(
                    f"Failed to load collection {genome_digest} from store '{name}' "
                    f"({self._rows[name].url}): {e}"
                ) from e
            try:
                store.pull_fhr(genome_digest)
            except Exception as e:  # noqa: BLE001 - a genome with no FHR is normal
                logger.debug(f"pull_fhr({genome_digest}) failed on store '{name}': {e}")
            self._loaded.add((name, genome_digest))
            return True

    def store_for_genome(
        self, genome_digest: GenomeDigest, store_name: str | None = None
    ) -> RefgetStore:
        """Return the store that owns ``genome_digest``, with it fully loaded.

        If ``store_name`` names a known store, it wins outright. Otherwise probe
        each store in priority order (stub-level, no collection fetch) and take
        the first that holds the digest. For a remote store, the genome's full
        collection and FHR sidecar are loaded on first use, so every caller
        gets a store it can read the genome from. Raises
        :class:`MissingGenomeError` if no store holds it.
        """
        if store_name and store_name in self._stores:
            if not self._ensure_loaded(store_name, genome_digest):
                raise MissingGenomeError(genome_digest)
            return self._stores[store_name]
        for name in self._order:
            try:
                held = self._holds(name, genome_digest)
            except Exception as e:  # noqa: BLE001 - a probe failure is not fatal
                logger.debug(f"Ownership probe failed on store '{name}': {e}")
                continue
            if held:
                self._ensure_loaded(name, genome_digest)
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
