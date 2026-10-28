"""Unit tests for the federated RefgetStoreRouter."""

import threading
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from refgenie.core.store_router import RefgetStoreRouter, make_store
from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError

pytestmark = pytest.mark.unit


class _FakeStore:
    """A RefgetStore stand-in holding collection stubs, like ``open_remote``.

    ``get_collection_metadata`` answers for any held digest (stub-level);
    ``is_collection_loaded`` is True only after ``load_collection``.
    """

    def __init__(self, digests, *, fhr_fails=False, load_delay=0.0):
        self._digests = set(digests)
        self._loaded = set()
        self._fhr = set()
        self._fhr_fails = fhr_fails
        self._load_delay = load_delay
        self.load_calls: list[str] = []
        self.fhr_calls: list[str] = []
        self.metadata_calls: list[str] = []

    def get_collection_metadata(self, digest):
        self.metadata_calls.append(digest)
        return {"digest": digest} if digest in self._digests else None

    def is_collection_loaded(self, digest):
        return digest in self._loaded

    def load_collection(self, digest):
        self.load_calls.append(digest)
        if digest not in self._digests:
            raise OSError("HTTP 403")
        time.sleep(self._load_delay)
        self._loaded.add(digest)

    def pull_fhr(self, digest=None):
        self.fhr_calls.append(digest)
        if self._fhr_fails:
            raise OSError("no FHR sidecar")
        self._fhr.add(digest)

    def get_collection_level2(self, digest):
        assert digest in self._loaded, "level 2 read from an unloaded stub"
        return {"names": ["chr1"], "lengths": [1], "sequences": ["SQ.x"]}

    def get_fhr_metadata(self, digest):
        return {"fhr": digest} if digest in self._fhr else None


class _FakeDiskStore(_FakeStore):
    """An on-disk store: everything is loaded already."""

    def __init__(self, digests):
        super().__init__(digests)
        self._loaded = set(digests)

    def get_collection_metadata(self, digest):  # pragma: no cover - must not be called
        raise AssertionError("on_disk stores are probed by is_collection_loaded")


def _router(monkeypatch, tmp_path, stores_by_url, rows):
    monkeypatch.setattr(
        "refgenie.core.store_router.make_store", lambda store, cache_dir: stores_by_url[store.url]
    )
    return RefgetStoreRouter(rows, tmp_path)


@pytest.fixture
def stores():
    return {
        "url://jungle": _FakeStore({"digA", "shared"}),
        "url://vgp": _FakeStore({"digB", "shared"}),
    }


@pytest.fixture
def router(monkeypatch, tmp_path, stores):
    """A router over two fake remote stores: 'jungle' (prio 10) and 'vgp' (prio 20)."""
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
    """An explicit, known store_name wins outright -- no probing."""
    # 'shared' is in both; explicit vgp must return vgp even though
    # jungle has higher priority.
    assert router.store_for_genome("shared", store_name="vgp") is router.get_store("vgp")


def test_store_for_genome_probe_fallback(router):
    """With no store_name, probe by digest in priority order."""
    assert router.store_for_genome("digB") is router.get_store("vgp")
    assert router.store_for_genome("digA") is router.get_store("jungle")
    # 'shared' is in both; the highest-priority store wins the probe.
    assert router.store_for_genome("shared") is router.get_store("jungle")


def test_store_for_genome_probes_stubs_by_metadata(router, stores):
    """Ownership is probed via get_collection_metadata, which works on stubs."""
    router.store_for_genome("digB")
    assert stores["url://jungle"].metadata_calls[0] == "digB"
    assert stores["url://jungle"].load_calls == []
    assert stores["url://vgp"].load_calls == ["digB"]


def test_store_for_genome_unknown_name_falls_back_to_probe(router):
    """An unknown store_name is ignored; the probe still resolves."""
    assert router.store_for_genome("digA", store_name="nope") is router.get_store("jungle")


def test_store_for_genome_missing_raises_without_loading(router, stores):
    """A digest no store holds raises MissingGenomeError before any fetch."""
    with pytest.raises(MissingGenomeError):
        router.store_for_genome("does-not-exist")
    with pytest.raises(MissingGenomeError):
        router.store_for_genome("does-not-exist", store_name="vgp")
    assert stores["url://jungle"].load_calls == []
    assert stores["url://vgp"].load_calls == []


def test_metadata_reads_return_none_when_unowned(router):
    """get_collection_level2/get_fhr_metadata return None for an unowned digest."""
    assert router.get_collection_level2("does-not-exist") is None
    assert router.get_fhr_metadata("does-not-exist") is None


def test_pass_throughs_load_lazily_once(router, stores):
    """The first routed read loads collection + FHR; later reads reuse them."""
    vgp = stores["url://vgp"]
    assert router.get_fhr_metadata("digB", store_name="vgp") == {"fhr": "digB"}
    assert router.get_collection_level2("digB")["names"] == ["chr1"]
    router.store_for_genome("digB")
    assert vgp.load_calls == ["digB"]
    assert vgp.fhr_calls == ["digB"]


def test_concurrent_first_requests_load_once(monkeypatch, tmp_path):
    """Eight threads asking for one genome at once trigger exactly one load."""
    store = _FakeStore({"digA"}, load_delay=0.05)
    rows = [Store(name="jungle", url="url://jungle", type=StoreType.remote, priority=10)]
    router = _router(monkeypatch, tmp_path, {"url://jungle": store}, rows)
    start = threading.Barrier(8)

    def hit(_):
        start.wait()
        return router.get_collection_level2("digA")

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(hit, range(8)))
    assert all(r is not None for r in results)
    assert store.load_calls == ["digA"]
    assert store.fhr_calls == ["digA"]


def test_pull_fhr_failure_still_marks_loaded(monkeypatch, tmp_path):
    """A genome without an FHR sidecar is loaded once and serves level 2."""
    store = _FakeStore({"digA"}, fhr_fails=True)
    rows = [Store(name="jungle", url="url://jungle", type=StoreType.remote, priority=10)]
    router = _router(monkeypatch, tmp_path, {"url://jungle": store}, rows)
    assert router.get_fhr_metadata("digA") is None
    assert router.get_collection_level2("digA") is not None
    assert store.load_calls == ["digA"]
    assert store.fhr_calls == ["digA"]


def test_load_failure_propagates(monkeypatch, tmp_path):
    """A network error while loading a held collection is an error, not a miss."""
    store = _FakeStore({"digA"})
    store.load_collection = lambda digest: (_ for _ in ()).throw(OSError("timeout"))
    rows = [Store(name="jungle", url="url://jungle", type=StoreType.remote, priority=10)]
    router = _router(monkeypatch, tmp_path, {"url://jungle": store}, rows)
    with pytest.raises(RuntimeError, match="Failed to load collection digA"):
        router.get_collection_level2("digA")


def test_on_disk_store_skips_lazy_path(monkeypatch, tmp_path):
    """on_disk stores are probed by is_collection_loaded and never lazily loaded."""
    disk = _FakeDiskStore({"digL"})
    rows = [Store(name="local", url="/tmp/x", type=StoreType.on_disk, priority=0)]
    router = _router(monkeypatch, tmp_path, {"/tmp/x": disk}, rows)
    assert router.store_for_genome("digL") is disk
    assert router.store_for_genome("digL", store_name="local") is disk
    assert disk.load_calls == []
    assert disk.fhr_calls == []


def test_make_store_remote_skips_eager_loads(monkeypatch, tmp_path):
    """Opening a remote store fetches the collection index (stubs) only."""
    calls = []

    class _Opened:
        def set_quiet(self, quiet):
            calls.append("set_quiet")

        def pull_aliases(self):  # pragma: no cover - must not be called
            calls.append("pull_aliases")

        def load_all_collections(self):  # pragma: no cover - must not be called
            calls.append("load_all_collections")

        def pull_fhr(self, *a):  # pragma: no cover - must not be called
            calls.append("pull_fhr")

    monkeypatch.setattr(
        "refgenie.core.store_router.RefgetStore.open_remote",
        lambda cache, url: calls.append("open_remote") or _Opened(),
    )
    row = Store(name="jungle", url="s3://bucket/jungle", type=StoreType.remote, priority=10)
    make_store(row, tmp_path / "jungle")
    assert calls == ["open_remote", "set_quiet"]
