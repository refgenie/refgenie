"""Unit tests for the federated RefgetStoreRouter."""

import pytest

from refgenie.core.store_router import RefgetStoreRouter
from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError


class _FakeStore:
    """A RefgetStore stand-in that knows which collection digests it holds."""

    def __init__(self, digests):
        self._digests = set(digests)

    # make_store calls these on open; no-ops here.
    def set_quiet(self, quiet):  # noqa: D401
        pass

    def load_all_collections(self):
        pass

    def pull_aliases(self, *a, **k):
        pass

    def pull_fhr(self, *a, **k):
        pass

    def is_collection_loaded(self, digest):
        return digest in self._digests


@pytest.fixture
def router(monkeypatch, tmp_path):
    """A router over two fake stores: 'jungle' (prio 10) and 'vgp' (prio 20)."""
    stores_by_url = {
        "url://jungle": _FakeStore({"digA", "shared"}),
        "url://vgp": _FakeStore({"digB", "shared"}),
    }

    def fake_make_store(store, cache_dir):
        return stores_by_url[store.url]

    monkeypatch.setattr("refgenie.core.store_router.make_store", fake_make_store)
    rows = [
        Store(name="vgp", url="url://vgp", type=StoreType.remote, priority=20),
        Store(name="jungle", url="url://jungle", type=StoreType.remote, priority=10),
    ]
    return RefgetStoreRouter(rows, tmp_path)


def test_priority_order(router):
    """Stores are ordered by ascending priority integer."""
    assert router.names == ["jungle", "vgp"]


def test_default_store_is_highest_priority(router):
    """default_store is the lowest-priority-integer store."""
    assert router.default_store is router.get_store("jungle")


def test_store_for_genome_explicit_name(router):
    """An explicit, known store_name wins outright -- no probing."""
    # 'shared' is loaded in both; explicit vgp must return vgp even though
    # jungle has higher priority.
    assert router.store_for_genome("shared", store_name="vgp") is router.get_store("vgp")


def test_store_for_genome_probe_fallback(router):
    """With no store_name, probe by digest in priority order."""
    assert router.store_for_genome("digB") is router.get_store("vgp")
    assert router.store_for_genome("digA") is router.get_store("jungle")
    # 'shared' is in both; the highest-priority store wins the probe.
    assert router.store_for_genome("shared") is router.get_store("jungle")


def test_store_for_genome_unknown_name_falls_back_to_probe(router):
    """An unknown store_name is ignored; the probe still resolves."""
    assert router.store_for_genome("digA", store_name="nope") is router.get_store("jungle")


def test_store_for_genome_missing_raises(router):
    """A digest no store holds raises MissingGenomeError."""
    with pytest.raises(MissingGenomeError):
        router.store_for_genome("does-not-exist")


def test_metadata_reads_return_none_when_unowned(router, monkeypatch):
    """get_collection_level2/get_fhr_metadata return None for an unowned digest."""
    assert router.get_collection_level2("does-not-exist") is None
    assert router.get_fhr_metadata("does-not-exist") is None
