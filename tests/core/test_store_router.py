"""Unit tests for the federated RefgetStoreRouter."""

import pytest

from refgenie.core.store_router import RefgetStoreRouter, make_store
from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError

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
    """An explicit, known store_name wins outright."""
    assert router.store_for_genome("shared", store_name="vgp") is router.get_store("vgp")


def test_store_for_genome_probe_fallback(router):
    """With no store_name, the first store holding the digest wins."""
    assert router.store_for_genome("digB") is router.get_store("vgp")
    assert router.store_for_genome("digA") is router.get_store("jungle")
    assert router.store_for_genome("shared") is router.get_store("jungle")


def test_store_for_genome_unknown_name_falls_back_to_probe(router):
    """An unknown store_name is ignored; the probe still resolves."""
    assert router.store_for_genome("digA", store_name="nope") is router.get_store("jungle")


def test_store_for_genome_missing_raises(router):
    """A digest no store holds raises MissingGenomeError."""
    with pytest.raises(MissingGenomeError):
        router.store_for_genome("does-not-exist")


def test_metadata_reads(router):
    """Pass-throughs read from the owning store, or return None when unowned."""
    assert router.get_fhr_metadata("digB") == {"fhr": "digB"}
    assert router.get_collection_level2("digA")["names"] == ["chr1"]
    assert router.get_collection_level2("does-not-exist") is None
    assert router.get_fhr_metadata("does-not-exist") is None


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
