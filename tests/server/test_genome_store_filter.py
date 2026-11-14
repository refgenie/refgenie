"""``GET /v4/genomes?store=...``: filter by store membership, not ownership.

Builds real on-disk RefgetStores and syncs them, so this is ``component`` tier.
"""

import pytest
from sqlmodel import Session

from refgenie.db.tables import AssetGroup

from tests.support.factories import build_on_disk_store, create_asset_class, server_with_stores
from tests.support.server import make_server_client, needs_server

pytestmark = [pytest.mark.component, needs_server]


@pytest.fixture
def world(tmp_path, fixtures_path):
    """X in jungle only, Y in vgp only, Z in both (owned by jungle). X has assets."""
    jungle, vgp = tmp_path / "jungle", tmp_path / "vgp"
    x = build_on_disk_store(jungle, fixtures_path / "rCRSd.fa", "rCRSd")
    z = build_on_disk_store(jungle, fixtures_path / "t7.fa", "t7")
    y = build_on_disk_store(vgp, fixtures_path / "demo.fa", "demo")
    assert build_on_disk_store(vgp, fixtures_path / "t7.fa", "t7") == z
    rg = server_with_stores(("jungle", jungle, 10), ("vgp", vgp, 20))
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
