"""Unit tests for the federated RefgetStoreRouter."""

import pytest

from refgenie.core.store_router import RefgetStoreRouter, make_store
from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError, UnservedStoreError

pytestmark = pytest.mark.unit


class _FakeStore:
    """A fully loaded store stand-in, like the router holds in server mode."""

    def __init__(self, digests, fhr=()):
        self._digests = set(digests)
        self._fhr = set(fhr)

    def is_collection_loaded(self, digest):
        return digest in self._digests

    def get_collection_level2(self, digest):
        assert digest in self._digests
        return {"names": ["chr1"], "lengths": [1], "sequences": ["SQ.x"]}

    def get_fhr_metadata(self, digest):
        return {"fhr": digest} if digest in self._fhr else None

    def list_collections(self, page, page_size):
        ordered = sorted(self._digests)
        chunk = ordered[page * page_size : (page + 1) * page_size]
        return {"results": [{"digest": d} for d in chunk]}


def _router(monkeypatch, tmp_path, stores_by_url, rows, readonly=True):
    seen = []

    def fake_make_store(store, cache_dir, *, readonly):
        seen.append(readonly)
        return stores_by_url[store.url]

    monkeypatch.setattr("refgenie.core.store_router.make_store", fake_make_store)
    router = RefgetStoreRouter(rows, tmp_path, readonly=readonly)
    assert seen == [readonly] * len(rows)
    return router


@pytest.fixture
def stores():
    return {
        "url://jungle": _FakeStore({"digA", "shared"}, fhr={"digA"}),
        "url://vgp": _FakeStore({"digB", "shared"}, fhr={"digB"}),
    }


@pytest.fixture
def router(monkeypatch, tmp_path, stores):
    """A router over two fake stores: 'jungle' (prio 10) and 'vgp' (prio 20)."""
    rows = [
        Store(name="vgp", url="url://vgp", type=StoreType.remote, priority=20),
        Store(name="jungle", url="url://jungle", type=StoreType.remote, priority=10),
    ]
    return _router(monkeypatch, tmp_path, stores, rows)


def test_priority_order(router):
    """Stores are ordered by ascending priority integer."""
    assert router.names == ["jungle", "vgp"]


def test_default_store_is_highest_priority(router):
    """default_store is the lowest-priority-integer store."""
    assert router.default_store is router.get_store("jungle")


def test_store_for_genome_explicit_name(router):
    """The genome's store_name picks its store."""
    assert router.store_for_genome("shared", "vgp") is router.get_store("vgp")
    assert router.store_for_genome("shared", "jungle") is router.get_store("jungle")


def test_store_for_genome_never_probes(router):
    """Another store holding the digest is never used."""
    with pytest.raises(MissingGenomeError):
        router.store_for_genome("digB", "jungle")


def test_store_for_genome_unknown_store_raises(router):
    """An unserved store_name raises UnservedStoreError, naming the store."""
    with pytest.raises(UnservedStoreError, match="'plantref', which this server does not serve"):
        router.store_for_genome("digA", "plantref")


def test_store_for_genome_no_store_needs_local(router):
    """No store_name means the node's own store, which this router lacks."""
    with pytest.raises(UnservedStoreError, match="has no store"):
        router.store_for_genome("digA", None)


def test_store_for_genome_no_store_reads_local(monkeypatch, tmp_path, stores):
    """No store_name reads the node's own 'local' store when there is one."""
    stores["url://local"] = _FakeStore({"mine"})
    rows = [
        Store(name="local", url="url://local", type=StoreType.on_disk, priority=0),
        Store(name="vgp", url="url://vgp", type=StoreType.remote, priority=20),
    ]
    router = _router(monkeypatch, tmp_path, stores, rows)
    assert router.store_for_genome("mine", None) is router.get_store("local")
    assert router.served_store_names() == (["vgp"], True)


def test_sole_store_reads_every_genome(monkeypatch, tmp_path, stores):
    """Local mode: the one store holds every genome, whatever its store_name."""
    rows = [Store(name="local", url="url://jungle", type=StoreType.on_disk, priority=0)]
    monkeypatch.setattr(
        "refgenie.core.store_router.make_store", lambda store, cache, *, readonly: stores[store.url]
    )
    router = RefgetStoreRouter(rows, tmp_path, readonly=False, sole_store=True)
    assert router.store_for_genome("digA", "vgp") is router.get_store("local")
    assert router.store_for_genome("digA", None) is router.get_store("local")
    assert router.serves("anything") and router.served_store_names() is None


def test_served_store_names(router):
    assert router.served_store_names() == (["jungle", "vgp"], False)
    assert router.serves("vgp") and not router.serves("plantref") and not router.serves(None)


def test_metadata_reads(router):
    """Pass-throughs read from the genome's store."""
    assert router.get_fhr_metadata("digB", "vgp") == {"fhr": "digB"}
    assert router.get_collection_level2("digA", "jungle")["names"] == ["chr1"]
    with pytest.raises(MissingGenomeError):
        router.get_collection_level2("does-not-exist", "jungle")
    assert router.get_fhr_metadata("does-not-exist", "jungle") is None


def test_readonly_flag_reaches_make_store(monkeypatch, tmp_path, stores):
    """The router passes its readonly flag to every make_store call."""
    rows = [Store(name="jungle", url="url://jungle", type=StoreType.remote, priority=10)]
    _router(monkeypatch, tmp_path, stores, rows, readonly=False)


class _Opened:
    """What RefgetStore.open_remote/on_disk returns, recording calls."""

    def __init__(self, calls):
        self.calls = calls
        self.frozen = object()

    def set_quiet(self, quiet):
        self.calls.append("set_quiet")

    def load_all_collections(self):
        self.calls.append("load_all_collections")

    def pull_fhr(self, *a):
        self.calls.append("pull_fhr")
        return {"pulled": 2, "skipped": 0, "not_found": 1, "conflicts": 0}

    def pull_aliases(self):  # pragma: no cover - must not be called
        self.calls.append("pull_aliases")

    def into_readonly(self):
        self.calls.append("into_readonly")
        return _Frozen()


class _Frozen:
    def stats(self):
        return {"n_collections": "3"}


def _patch_open(monkeypatch, calls, opened):
    monkeypatch.setattr(
        "refgenie.core.store_router.RefgetStore.open_remote",
        lambda cache, url: calls.append("open_remote") or opened,
    )
    monkeypatch.setattr(
        "refgenie.core.store_router.RefgetStore.on_disk",
        lambda path: calls.append("on_disk") or opened,
    )


def test_make_store_remote_readonly_loads_everything(monkeypatch, tmp_path):
    """Server mode: a remote store loads every collection and FHR once, frozen."""
    calls = []
    opened = _Opened(calls)
    _patch_open(monkeypatch, calls, opened)
    row = Store(name="jungle", url="s3://bucket/jungle", type=StoreType.remote, priority=10)
    result = make_store(row, tmp_path / "jungle", readonly=True)
    assert isinstance(result, _Frozen)
    assert calls == [
        "open_remote",
        "set_quiet",
        "load_all_collections",
        "pull_fhr",
        "into_readonly",
    ]


def test_make_store_remote_readonly_failure_raises(monkeypatch, tmp_path):
    """A store that cannot load fails loudly instead of serving half a catalog."""
    calls = []
    opened = _Opened(calls)
    opened.load_all_collections = lambda: (_ for _ in ()).throw(OSError("HTTP 403"))
    _patch_open(monkeypatch, calls, opened)
    row = Store(name="jungle", url="s3://bucket/jungle", type=StoreType.remote, priority=10)
    with pytest.raises(OSError, match="403"):
        make_store(row, tmp_path / "jungle", readonly=True)


def test_make_store_on_disk_readonly_freezes(monkeypatch, tmp_path):
    """Server mode: the node's own on-disk store is loaded and frozen too."""
    calls = []
    _patch_open(monkeypatch, calls, _Opened(calls))
    row = Store(name="local", url=str(tmp_path / "s"), type=StoreType.on_disk, priority=0)
    assert isinstance(make_store(row, tmp_path / "c", readonly=True), _Frozen)
    assert calls == ["on_disk", "set_quiet", "load_all_collections", "into_readonly"]


def test_make_store_editable(monkeypatch, tmp_path):
    """readonly=False: on_disk stays editable; remote opens stubs only."""
    calls = []
    opened = _Opened(calls)
    _patch_open(monkeypatch, calls, opened)
    disk = Store(name="local", url=str(tmp_path / "s"), type=StoreType.on_disk, priority=0)
    assert make_store(disk, tmp_path / "c", readonly=False) is opened
    assert calls == ["on_disk", "set_quiet", "load_all_collections"]
    calls.clear()
    remote = Store(name="jungle", url="s3://bucket/jungle", type=StoreType.remote, priority=10)
    assert make_store(remote, tmp_path / "jungle", readonly=False) is opened
    assert calls == ["open_remote", "set_quiet"]


def test_collection_digests_is_membership(router):
    """Each store's set holds every digest it lists, shared ones included."""
    assert router.collection_digests("jungle") == {"digA", "shared"}
    assert router.collection_digests("vgp") == {"digB", "shared"}


def test_collection_digests_pages_past_one_page(monkeypatch, tmp_path):
    digests = {f"d{i:04d}" for i in range(2500)}
    rows = [Store(name="big", url="url://big", type=StoreType.remote, priority=1)]
    router = _router(monkeypatch, tmp_path, {"url://big": _FakeStore(digests)}, rows)
    assert router.collection_digests("big") == digests


def test_collection_digests_cached_when_readonly(router):
    assert router.collection_digests("vgp") is router.collection_digests("vgp")


def test_collection_digests_recomputed_when_editable(monkeypatch, tmp_path, stores):
    rows = [Store(name="vgp", url="url://vgp", type=StoreType.remote, priority=20)]
    router = _router(monkeypatch, tmp_path, stores, rows, readonly=False)
    first = router.collection_digests("vgp")
    stores["url://vgp"]._digests.add("new")
    assert router.collection_digests("vgp") == first | {"new"}


def test_collection_digests_unknown_store_raises(router):
    with pytest.raises(KeyError):
        router.collection_digests("nope")


@pytest.mark.component
def test_collection_digests_real_on_disk_stores(tmp_path, fixtures_path):
    """Against frozen on-disk stores, the set is exactly what was ingested."""
    from refget.store import RefgetStore

    def build(path, fastas):
        store = RefgetStore.on_disk(str(path))
        return {store.add_sequence_collection_from_fasta(str(f))[0].digest for f in fastas}

    a = build(tmp_path / "a", [fixtures_path / "rCRSd.fa", fixtures_path / "t7.fa"])
    b = build(tmp_path / "b", [fixtures_path / "demo.fa"])
    rows = [
        Store(name="a", url=str(tmp_path / "a"), type=StoreType.on_disk, priority=1),
        Store(name="b", url=str(tmp_path / "b"), type=StoreType.on_disk, priority=2),
    ]
    router = RefgetStoreRouter(rows, tmp_path / "cache", readonly=True)
    assert router.collection_digests("a") == a
    assert router.collection_digests("b") == b


@pytest.mark.component
@pytest.mark.parametrize("scan", [False, True])
def test_collection_sequences_bulk_and_cached(tmp_path, fixtures_path, monkeypatch, scan):
    """Rows match level 2 with metadata and aliases, by either alias path, and are cached."""
    from refget.store import RefgetStore

    store = RefgetStore.on_disk(str(tmp_path / "s"))
    digest = store.add_sequence_collection_from_fasta(str(fixtures_path / "demo.fa"))[0].digest
    level2 = store.get_collection_level2(digest)
    first = level2["sequences"][0].removeprefix("SQ.")
    store.add_sequence_alias("insdc", "XY000001.1", first)
    rows = [Store(name="s", url=str(tmp_path / "s"), type=StoreType.on_disk, priority=1)]
    router = RefgetStoreRouter(rows, tmp_path / "cache", readonly=True)
    monkeypatch.setattr(RefgetStoreRouter, "ALIAS_SCAN_MIN_SEQUENCES", 0 if scan else 10**9)

    result = router.get_collection_sequences(digest, "s")

    assert [r.name for r in result.sequences] == level2["names"]
    assert [r.length for r in result.sequences] == level2["lengths"]
    assert [r.sha512t24u for r in result.sequences] == [
        d.removeprefix("SQ.") for d in level2["sequences"]
    ]
    assert all(r.md5 and r.alphabet for r in result.sequences)
    assert [(a.namespace, a.alias) for a in result.sequences[0].aliases] == [
        ("insdc", "XY000001.1")
    ]
    assert all(not r.aliases for r in result.sequences[1:])
    assert result.store == "s"
    assert router.get_collection_sequences(digest, "s") is result
