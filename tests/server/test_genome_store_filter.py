"""``GET /v4/genomes?store=...``: filter by store membership, not ownership.

Builds real on-disk RefgetStores and syncs them, so this is ``component`` tier.
"""

import pytest
from sqlmodel import Session

from refgenie import Refgenie
from refgenie.db.tables import AssetGroup, StoreType

from tests.support.factories import create_asset_class, make_engine
from tests.support.server import make_server_client, needs_server

pytestmark = [pytest.mark.component, needs_server]


def _store(path, fastas):
    """An on-disk RefgetStore holding these FASTAs; returns their digests."""
    from refget.store import RefgetStore

    store = RefgetStore.on_disk(str(path))
    return [store.add_sequence_collection_from_fasta(str(f))[0].digest for f in fastas]


@pytest.fixture
def world(tmp_path, fixtures_path):
    """X in jungle only, Y in vgp only, Z in both (owned by jungle). X has assets."""
    from refgenie.cli.commands.store import _sync_one_store

    x, z = _store(tmp_path / "jungle", [fixtures_path / "rCRSd.fa", fixtures_path / "t7.fa"])
    y, z2 = _store(tmp_path / "vgp", [fixtures_path / "demo.fa", fixtures_path / "t7.fa"])
    assert z == z2
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True, server_mode=True)
    rg.database.create_tables()
    rg.store.add("jungle", str(tmp_path / "jungle"), StoreType.on_disk, priority=10)
    rg.store.add("vgp", str(tmp_path / "vgp"), StoreType.on_disk, priority=20)
    for srow in rg.store.enabled_stores():
        _sync_one_store(rg, srow, {"jungle": 10, "vgp": 20}, 100)
    with Session(rg.database_engine) as session:
        ac = create_asset_class(session, "fasta", ["files"])
        session.add(AssetGroup(name="fasta", description="", genome_digest=x, asset_class_id=ac.id))
        session.commit()
    assert rg.genome.get(z).store_name == "jungle"
    return rg, {"X": x, "Y": y, "Z": z}


def _digests(client, query):
    response = client.get(f"/v4/genomes?{query}")
    assert response.status_code == 200, response.text
    return {g["digest"] for g in response.json()["results"]}


def test_store_filter_is_membership(world):
    rg, d = world
    with make_server_client(rg) as client:
        # Z is owned by jungle but still present in vgp.
        assert _digests(client, "store=vgp") == {d["Y"], d["Z"]}
        assert _digests(client, "store=jungle") == {d["X"], d["Z"]}
        assert _digests(client, "store=vgp,jungle") == {d["X"], d["Y"], d["Z"]}
        assert _digests(client, "store=%20vgp%20,,") == {d["Y"], d["Z"]}


def test_store_filter_combines_with_has_assets(world):
    rg, d = world
    with make_server_client(rg) as client:
        assert _digests(client, "store=jungle&has_assets=true") == {d["X"]}
        assert _digests(client, "store=vgp&has_assets=true") == set()
        assert _digests(client, "store=vgp&has_assets=false") == {d["Y"], d["Z"]}


def test_unknown_store_is_422(world):
    rg, _ = world
    with make_server_client(rg) as client:
        response = client.get("/v4/genomes?store=nope")
    assert response.status_code == 422
    detail = response.json()["detail"]
    assert "nope" in detail and "jungle" in detail and "vgp" in detail


def test_listing_carries_assembly_level(world):
    rg, _ = world
    with make_server_client(rg) as client:
        results = client.get("/v4/genomes").json()["results"]
    assert results and all("assembly_level" in g for g in results)
