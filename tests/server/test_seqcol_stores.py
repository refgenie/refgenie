"""The server serves whatever stores the registry lists, and follows it at runtime.

A catalog import writes the `store` table; these tests stand in for one with
``rg.store.add`` / ``apply_declaration`` and then call the refresh the import
runs (``app.state.seqcol.refresh()``). The app is built once per test: a store
change must not need a new app.

They build real on-disk RefgetStores, so they are ``component`` tier.
"""

import pytest

from refgenie.db.tables import StoreType

from tests.support.factories import build_on_disk_store, make_engine
from tests.support.server import make_server_app, needs_server

pytestmark = [pytest.mark.component, needs_server]

FHR = {
    "genome": "rCRSd",
    "taxon": {"name": "Homo sapiens", "uri": "https://identifiers.org/taxonomy:9606"},
}


def _server_rg():
    from refgenie import Refgenie

    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True, server_mode=True)
    rg.database.create_tables()
    return rg


def _client(rg):
    from fastapi.testclient import TestClient

    return TestClient(make_server_app(rg))


class _Stores:
    """Three on-disk stores, one genome each; registered on demand."""

    PRIORITY = {"jungle": 10, "vgp": 20, "plantref": 30}

    def __init__(self, tmp_path, fixtures_path):
        from refget.store import RefgetStore

        from refgenie.managers.genome import GenomeManager

        self.paths = {name: tmp_path / name for name in self.PRIORITY}
        self.digest = {
            "jungle": build_on_disk_store(
                self.paths["jungle"], fixtures_path / "rCRSd.fa", "rCRSd"
            ),
            "vgp": build_on_disk_store(self.paths["vgp"], fixtures_path / "t7.fa", "t7"),
            "plantref": build_on_disk_store(
                self.paths["plantref"], fixtures_path / "demo.fa", "demo"
            ),
        }
        RefgetStore.on_disk(str(self.paths["jungle"])).set_fhr_metadata(
            self.digest["jungle"], GenomeManager._fhr_metadata_from_dict(FHR)
        )

    def row(self, name, enabled=True):
        """A store declaration, as a catalog import would write it."""
        return {
            "name": name,
            "url": str(self.paths[name]),
            "type": StoreType.on_disk.value,
            "priority": self.PRIORITY[name],
            "enabled": enabled,
        }

    def register(self, rg, *names):
        """Add the stores and their genomes, as a catalog import does."""
        from refgenie.cli.commands.store import _sync_one_store

        rg.store.apply_declaration([self.row(n) for n in names])
        for name in names:
            _sync_one_store(rg, rg.store.get(name), self.PRIORITY, 100)


@pytest.fixture
def stores(tmp_path, fixtures_path):
    return _Stores(tmp_path, fixtures_path)


def _served_stores(client):
    """Store names in /service-info's seqcol.services, default first."""
    services = client.get("/service-info").json()["seqcol"]["services"]
    return [s.get("store") for s in services]


def _serves_collection(client, prefix, digest):
    return client.get(f"{prefix}/collection/{digest}?level=2").status_code == 200


class TestNoStores:
    def test_db_mode_seqcol(self):
        with _client(_server_rg()) as client:
            assert client.app.state.seqcol.mounts == ()
            info = client.get("/seqcol/service-info").json()
            assert info["id"] == "org.refgenie.seqcol"
            assert info["seqcol"]["refget_store"] == {"enabled": False}
            assert _served_stores(client) == [None]
            # The DB-backed router answers; the collection simply is not there.
            assert client.get("/seqcol/list/collection").status_code == 200


class TestStoresAreServed:
    def test_each_store_serves_seqcol_and_v4(self, stores):
        rg = _server_rg()
        stores.register(rg, "jungle", "vgp")
        with _client(rg) as client:
            assert _served_stores(client) == ["jungle", "vgp"]
            for name, prefix in (("jungle", "/seqcol"), ("vgp", "/seqcol/vgp")):
                digest = stores.digest[name]
                assert _serves_collection(client, prefix, digest)
                info = client.get(f"{prefix}/service-info").json()
                assert info["seqcol"]["refget_store"]["url"] == str(stores.paths[name])
                assert client.get(f"/v4/genomes/{digest}").status_code == 200
            alias = client.get("/v4/aliases/rCRSd").json()
            assert alias["collection"]["names"]
            assert alias["fhr"]["genome"] == "rCRSd"
            assert client.get(f"/v4/genomes/{stores.digest['jungle']}").json()["fhr"]
            assert client.get("/v4/aliases/t7").status_code == 200
            # The default store is not also served under its own name.
            assert client.get("/seqcol/jungle/service-info").status_code == 404


class TestStoreSetChangesAtRuntime:
    def test_added_store_is_served_after_refresh(self, stores):
        rg = _server_rg()
        stores.register(rg, "jungle")
        with _client(rg) as client:
            app = client.app
            plant = stores.digest["plantref"]
            assert not _serves_collection(client, "/seqcol/plantref", plant)
            assert client.get("/v4/aliases/demo").status_code == 404

            stores.register(rg, "plantref")
            assert app.state.seqcol.refresh() is True

            assert client.app is app
            assert _served_stores(client) == ["jungle", "plantref"]
            assert _serves_collection(client, "/seqcol/plantref", plant)
            assert _serves_collection(client, "/seqcol", stores.digest["jungle"])
            assert client.get("/v4/aliases/demo").status_code == 200
            assert rg.store_router.names == ["jungle", "plantref"]

    def test_disabled_store_stops_being_served(self, stores):
        rg = _server_rg()
        stores.register(rg, "jungle", "vgp")
        with _client(rg) as client:
            phage = stores.digest["vgp"]
            assert _serves_collection(client, "/seqcol/vgp", phage)

            rg.store.apply_declaration([stores.row("vgp", enabled=False)])
            assert client.app.state.seqcol.refresh() is True

            assert _served_stores(client) == ["jungle"]
            assert client.get("/seqcol/vgp/service-info").status_code == 404
            assert not _serves_collection(client, "/seqcol", phage)
            assert client.get("/v4/aliases/t7").status_code == 404
            assert rg.store_router.names == ["jungle"]

    def test_new_default_store_takes_the_default_mount(self, stores):
        rg = _server_rg()
        stores.register(rg, "vgp")
        with _client(rg) as client:
            assert _serves_collection(client, "/seqcol", stores.digest["vgp"])
            stores.register(rg, "jungle")  # priority 10 beats vgp's 20
            client.app.state.seqcol.refresh()
            assert _served_stores(client) == ["jungle", "vgp"]
            assert _serves_collection(client, "/seqcol", stores.digest["jungle"])
            assert _serves_collection(client, "/seqcol/vgp", stores.digest["vgp"])

    def test_disabling_every_store_falls_back_to_db_mode(self, stores):
        rg = _server_rg()
        stores.register(rg, "jungle")
        with _client(rg) as client:
            rg.store.apply_declaration([stores.row("jungle", enabled=False)])
            assert client.app.state.seqcol.refresh() is True
            assert client.app.state.seqcol.mounts == ()
            info = client.get("/seqcol/service-info").json()
            assert info["seqcol"]["refget_store"] == {"enabled": False}

    def test_unchanged_store_set_keeps_everything(self, stores):
        rg = _server_rg()
        stores.register(rg, "jungle")
        with _client(rg) as client:
            router = rg.store_router
            mounts = client.app.state.seqcol.mounts
            assert client.app.state.seqcol.refresh() is False
            assert rg.store_router is router
            assert client.app.state.seqcol.mounts is mounts

    def test_a_store_that_fails_to_load_keeps_the_old_set(self, stores):
        rg = _server_rg()
        stores.register(rg, "jungle")
        with _client(rg) as client:
            router = rg.store_router
            # Nothing listens on port 9, so opening this store fails.
            rg.store.add("broken", "http://127.0.0.1:9/", StoreType.remote, priority=50)
            with pytest.raises(Exception):
                client.app.state.seqcol.refresh()
            assert rg.store_router is router
            assert _served_stores(client) == ["jungle"]
            assert _serves_collection(client, "/seqcol", stores.digest["jungle"])

    def test_an_unreadable_registry_keeps_the_current_stores(self, stores, monkeypatch):
        rg = _server_rg()
        stores.register(rg, "jungle")
        with _client(rg) as client:

            def _fail():
                raise RuntimeError("database went away")

            monkeypatch.setattr(rg.store, "enabled_stores", _fail)
            assert client.app.state.seqcol.refresh() is False
            assert _served_stores(client) == ["jungle"]
            assert _serves_collection(client, "/seqcol", stores.digest["jungle"])
