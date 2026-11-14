"""Federated store router.

Opens one refget store per enabled :class:`~refgenie.db.tables.Store` row, in
priority order, and dispatches reads to the store that owns a genome. In server
mode (``readonly=True``) each store is loaded in full and frozen, so no request
touches the network. A change to the set of stores builds a new router
(``Refgenie.reload_store_router``); a change to a store's contents needs a
restart. In local mode
(``readonly=False``) it holds the editable on-disk store.
"""

import threading
import time
from collections import OrderedDict
from pathlib import Path

from refget.store import ReadonlyRefgetStore, RefgetStore

from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError, UnservedStoreError
from refgenie.logger import logger
from refgenie.managers.store import LOCAL_STORE_NAME
from refgenie.models import AliasEntry, CollectionSequences, GenomeDigest, SequenceRow

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


class _LevelTwoRecord:
    """Sequence metadata from a level-2 row plus the store's sequence record."""

    def __init__(self, name: str, length: int, sha512t24u: str, store):
        meta = store.get_sequence_metadata(sha512t24u)
        self.name = name
        self.length = length
        self.sha512t24u = sha512t24u
        self.md5 = meta.md5 if meta else None
        self.alphabet = meta.alphabet if meta else None
        self.description = meta.description if meta else None


class RefgetStoreRouter:
    """Holds one opened store per enabled :class:`Store` row, priority-ordered."""

    def __init__(
        self, stores: list[Store], cache_dir: Path, *, readonly: bool, sole_store: bool = False
    ):
        """
        Args:
            stores: The store rows to open.
            cache_dir: Where remote stores cache what they fetch.
            readonly: Load each store in full and freeze it (server mode).
            sole_store: The one store (there must be exactly one) holds every
                genome this node reads, so a genome's ``store_name`` does not
                pick it (local mode).
        """
        if sole_store and len(stores) != 1:
            raise ValueError("A sole-store router needs exactly one store")
        self._cache_dir = Path(cache_dir)
        self._readonly = readonly
        self._sole_store = sole_store
        #: (store, digest) -> rows, for frozen stores only (collections never change).
        self._sequences_cache: OrderedDict[tuple[str, str], CollectionSequences] = OrderedDict()
        self._sequences_lock = threading.Lock()

        #: Cached digest sets per store; filled only when the stores are frozen.
        self._digests: dict[str, frozenset[str]] = {}
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

    def collection_digests(self, name: str) -> frozenset[str]:
        """Digests of every collection in store ``name`` (membership, not ownership).

        A genome held by several stores is in each store's set, whichever store
        owns it. Results are cached when the router is readonly (frozen stores
        never change); editable stores are re-listed on each call. An unknown
        name raises :class:`KeyError`.
        """
        if name in self._digests:
            return self._digests[name]
        store = self.get_store(name)
        page_size = 1000
        digests: set[str] = set()
        page = 0
        while True:
            result = store.list_collections(page=page, page_size=page_size)
            items = (result or {}).get("results", [])
            for item in items:
                digest = (
                    item.get("digest") if isinstance(item, dict) else getattr(item, "digest", None)
                )
                if digest:
                    digests.add(digest)
            if len(items) < page_size:
                break
            page += 1
        frozen = frozenset(digests)
        if self._readonly:
            self._digests[name] = frozen
        return frozen

    # --- Dispatch -----------------------------------------------------------

    @property
    def sole_store(self) -> bool:
        """Whether one store holds every genome this router reads (local mode)."""
        return self._sole_store

    def serves(self, store_name: str | None) -> bool:
        """Whether a genome whose ``store_name`` is this is read through this router.

        ``None`` means a genome this node built itself, held by its own
        ``local`` store.
        """
        if self._sole_store:
            return True
        return (store_name or LOCAL_STORE_NAME) in self._stores

    def served_store_names(self) -> tuple[list[str], bool] | None:
        """The ``store_name`` values of the genomes this router serves.

        Returns ``(names, include_storeless)``, or ``None`` when every genome
        is served: in local mode, and on a server with no stores at all (DB
        seqcol), which routes nothing.
        """
        if self._sole_store or not self._order:
            return None
        names = [n for n in self._order if n != LOCAL_STORE_NAME]
        return names, LOCAL_STORE_NAME in self._stores

    def store_for_genome(
        self, genome_digest: GenomeDigest, store_name: str | None
    ) -> AnyRefgetStore:
        """Return the store that holds ``genome_digest``, with its collection loaded.

        ``store_name`` is the genome's ``store_name`` and the only thing that
        picks the store; no store is probed. ``None`` means the genome was
        built on this node and lives in its ``local`` store. In local mode the
        one on-disk store holds every genome, whatever its ``store_name``.

        Raises:
            UnservedStoreError: If this router does not serve that store.
            MissingGenomeError: If the store does not hold the collection.
        """
        if self._sole_store:
            store = self._stores[self._order[0]]
        else:
            store = self._stores.get(store_name or LOCAL_STORE_NAME)
            if store is None:
                raise UnservedStoreError(genome_digest, store_name)
        if not store.is_collection_loaded(genome_digest):
            raise MissingGenomeError(genome_digest)
        return store

    def _store_name_for(self, store: AnyRefgetStore) -> str:
        """The registered name of an opened store."""
        return next(name for name, s in self._stores.items() if s is store)

    # --- Genome-routed pass-throughs ---------------------------------------
    # Each takes the genome's ``store_name``, which picks the store.

    def get_sequence_by_name(self, genome_digest: GenomeDigest, name: str, store_name: str | None):
        """Fetch a sequence record from the store that holds ``genome_digest``."""
        return self.store_for_genome(genome_digest, store_name).get_sequence_by_name(
            genome_digest, name
        )

    def get_collection_level2(self, genome_digest: GenomeDigest, store_name: str | None):
        """Fetch level-2 collection data.

        Raises:
            MissingGenomeError: If the genome's store is not served or lacks it.
        """
        return self.store_for_genome(genome_digest, store_name).get_collection_level2(genome_digest)

    def get_fhr_metadata(self, genome_digest: GenomeDigest, store_name: str | None):
        """Fetch FHR metadata, or ``None`` if the genome's store lacks it or none is set."""
        try:
            store = self.store_for_genome(genome_digest, store_name)
        except MissingGenomeError:
            return None
        return store.get_fhr_metadata(genome_digest)

    #: Most collections :meth:`get_collection_sequences` keeps per router.
    SEQUENCES_CACHE_SIZE = 32

    def get_collection_sequences(
        self, genome_digest: GenomeDigest, store_name: str | None
    ) -> CollectionSequences:
        """A genome's sequences with metadata and aliases.

        Sequence records come from one ``get_collection`` call, and aliases
        from a per-store index built once, so the cost does not grow with a
        store call per sequence. Frozen stores never change, so results are
        kept in a small LRU cache keyed by (store, digest).

        Raises:
            MissingGenomeError: If the genome's store is not served or lacks it.
        """
        store = self.store_for_genome(genome_digest, store_name)
        name = self._store_name_for(store)
        key = (name, genome_digest)
        if self._readonly:
            with self._sequences_lock:
                cached = self._sequences_cache.get(key)
                if cached is not None:
                    self._sequences_cache.move_to_end(key)
                    return cached
        result = self._read_collection_sequences(name, store, genome_digest)
        if self._readonly:
            with self._sequences_lock:
                self._sequences_cache[key] = result
                self._sequences_cache.move_to_end(key)
                while len(self._sequences_cache) > self.SEQUENCES_CACHE_SIZE:
                    self._sequences_cache.popitem(last=False)
        return result

    #: From this many sequences up, aliases come from one pass over the
    #: store's alias tables rather than a store call per sequence. Each such
    #: call scans the tables (about 4 ms on the vgp store); one pass costs
    #: about 0.5 s however many sequences there are.
    ALIAS_SCAN_MIN_SEQUENCES = 100

    def _sequence_aliases(self, store: AnyRefgetStore, digests: set[str]):
        """``{bare sequence digest: aliases}`` for ``digests``, in store ``store``.

        Nothing is kept: a whole store's alias index would cost far more
        memory than the per-collection result cache.
        """
        if len(digests) < self.ALIAS_SCAN_MIN_SEQUENCES or not self._readonly:
            return {
                d: [AliasEntry(ns, a) for ns, a in store.get_aliases_for_sequence(d)]
                for d in digests
            }
        found: dict[str, list[AliasEntry]] = {}
        for namespace in store.list_sequence_alias_namespaces():
            for alias in store.list_sequence_aliases(namespace) or []:
                meta = store.get_sequence_metadata_by_alias(namespace, alias)
                if meta is not None and meta.sha512t24u in digests:
                    found.setdefault(meta.sha512t24u, []).append(AliasEntry(namespace, alias))
        return found

    def _read_collection_sequences(
        self, name: str, store: AnyRefgetStore, genome_digest: GenomeDigest
    ) -> CollectionSequences:
        """Build :meth:`get_collection_sequences`'s result from the store."""
        if self._readonly:
            # One call; a frozen store hands back records without sequence bytes.
            records = [record.metadata for record in store.get_collection(genome_digest).sequences]
        else:
            # An editable store's get_collection loads the bytes; read metadata only.
            level2 = store.get_collection_level2(genome_digest)
            records = [
                _LevelTwoRecord(name, length, digest.removeprefix("SQ."), store)
                for name, length, digest in zip(
                    level2["names"], level2["lengths"], level2["sequences"], strict=True
                )
            ]
        aliases = self._sequence_aliases(store, {meta.sha512t24u for meta in records})
        rows = []
        for meta in records:
            rows.append(
                SequenceRow(
                    name=meta.name,
                    length=meta.length,
                    sha512t24u=meta.sha512t24u,
                    md5=meta.md5,
                    alphabet=str(meta.alphabet) if meta.alphabet is not None else None,
                    description=meta.description,
                    aliases=list(aliases.get(meta.sha512t24u, ())),
                )
            )
        return CollectionSequences(
            digest=genome_digest,
            store=name,
            level1=dict(store.get_collection_level1(genome_digest)),
            sequences=rows,
            collection_aliases=[
                AliasEntry(ns, a) for ns, a in store.get_aliases_for_collection(genome_digest)
            ],
        )
