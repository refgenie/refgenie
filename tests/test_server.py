"""
Tests for the refgenie HTTP servers.

Three groups live here, because they are all "the server as an HTTP app":

* **v4 server endpoints** -- the alias and genomes endpoints, the three
  service-info documents, file serving, the /archives list and resolved
  serving_modes, over real on-disk builds.
* **core endpoints** -- the shared /v4 router over an empty in-memory DB with a
  stubbed refgenie. Pagination envelope/bounds and lookup 404s are parametrized
  across every paginated endpoint.
* **the data-channel router** -- /v4/data_channel/<channel>/... as a
  transparent, structure-preserving mirror of a channel's upstream index: the
  served index.yaml and the relative file paths are identical to what the
  upstream (e.g. the refgenie-registry GitHub Pages channel) exposes -- native
  layout, nested recipe paths, no filename rewriting or flattening.

Then the route-registration and middleware guards over both apps.

(RefgenieserverClient pagination lives with the sources/remote tests, not
here -- it exercises the client, not the app.)
"""

import re
from collections import Counter, defaultdict
from unittest.mock import MagicMock

import pytest
import yaml

from tests.helpers import (
    make_engine,
    make_local_client,
    make_server_app,
    make_server_client,
    make_server_rgc,
    register_fasta,
    requires_dash,
    requires_server,
    server_with_asset,
    stub_rgc,
)

# Server extras are optional; skip the whole module if any is missing.
requires_server()
requires_dash()

from fastapi.routing import APIRoute  # noqa: E402  (must follow the extras guard)

# The route-registration guards below walk a FastAPI app's *effective* routes:
# since 0.137, include_router no longer flattens child routes onto app.routes
# but hides them behind an opaque _IncludedRouter, so the prefixed
# path/methods/operation_id must be recovered with FastAPI's own walker.
# `iter_route_contexts` is the supported (non-underscore) walker get_openapi
# itself uses. It is a version-specific internal, so import it defensively: a
# FastAPI that does not expose it makes the guards skip (see
# `_effective_api_routes`) instead of crashing collection of this whole module
# -- and, via `from tests.test_server import ...`, of tests/test_web.py too.
try:
    from fastapi.routing import iter_route_contexts  # noqa: E402
except ImportError:  # pragma: no cover - depends on FastAPI internals
    iter_route_contexts = None
from fastapi.testclient import TestClient  # noqa: E402
from sqlmodel import Session, SQLModel  # noqa: E402

from refgenie import Refgenie  # noqa: E402
from refgenie.db.tables import Asset, AssetGroup, AssetName  # noqa: E402
from refgenie.managers.sources.manager import IndexFile, IndexFileSection  # noqa: E402
from refgenie.server.schemas import DataChannel  # noqa: E402
from refgenie.utils.build import get_dir_digest  # noqa: E402

# --- v4 server endpoints -----------------------------------------------------


@pytest.fixture
def rgc_with_alias(tmp_path, fixtures_path):
    """A catalog holding one genome (test_digest_12345) with two aliases."""
    return make_server_rgc(
        tmp_path, fixtures_path, genomes=[("test_digest_12345", ["rCRSd", "hg38_mito"])]
    )


@pytest.fixture
def client_with_alias(rgc_with_alias):
    """Client over ``rgc_with_alias``."""
    with make_server_client(rgc_with_alias) as client:
        yield client


@pytest.fixture
def client_no_alias(tmp_path, fixtures_path):
    """Client over a genome with no aliases."""
    rgc = make_server_rgc(tmp_path, fixtures_path, genomes=[("no_alias_digest_999", [])])
    with make_server_client(rgc) as client:
        yield client


@pytest.fixture
def client_plain(tmp_path, fixtures_path):
    """Client over an empty catalog (no genomes)."""
    with make_server_client(make_server_rgc(tmp_path, fixtures_path)) as client:
        yield client


@pytest.fixture
def server_app(tmp_path, fixtures_path):
    """The server app object (for route inspection)."""
    return make_server_app(make_server_rgc(tmp_path, fixtures_path))


@pytest.fixture(scope="module")
def service_info_client():
    """The module-global server app (DB-backed, no store URL) for service-info tests."""
    from refgenie.server.main import app

    with TestClient(app) as client:
        yield client


class TestAliasEndpoints:
    """GET /v4/aliases and GET /v4/aliases/{name}."""

    def test_list_aliases(self, client_with_alias):
        client = client_with_alias
        data = client.get("/v4/aliases").json()
        names = [item["name"] for item in data["items"]]
        assert "rCRSd" in names
        assert "hg38_mito" in names

    def test_list_aliases_filter_by_genome(self, client_with_alias):
        client = client_with_alias
        assert len(client.get("/v4/aliases?genome_digest=test_digest_12345").json()["items"]) == 2
        assert len(client.get("/v4/aliases?genome_digest=nonexistent").json()["items"]) == 0

    def test_list_aliases_is_paginated(self, client_with_alias):
        """The client's page loop needs `pagination` and honored offset/limit.

        An unpaginated list makes ServerSource.get_paginated repeat forever once a
        server has more aliases than one page.
        """
        client = client_with_alias
        assert client.get("/v4/aliases").json()["pagination"] == {
            "offset": 0,
            "limit": 100,
            "total": 2,
        }
        page = client.get("/v4/aliases?offset=1&limit=1").json()
        assert page["pagination"] == {"offset": 1, "limit": 1, "total": 2}
        assert len(page["items"]) == 1

    @pytest.mark.parametrize(
        "path",
        [
            "/v4/aliases?offset=-1",
            "/v4/aliases?limit=0",
            "/v4/aliases?limit=1001",
        ],
    )
    def test_list_aliases_rejects_bad_pagination(self, client_with_alias, path):
        client = client_with_alias
        assert client.get(path).status_code == 422

    def test_list_aliases_filter_by_name_and_search(self, client_with_alias):
        client = client_with_alias
        assert [i["name"] for i in client.get("/v4/aliases?name=rCRSd").json()["items"]] == [
            "rCRSd"
        ]
        found = client.get("/v4/aliases?q=mito&search_fields=name").json()
        assert [i["name"] for i in found["items"]] == ["hg38_mito"]
        missing = client.get("/v4/aliases?q=zzzznomatch").json()
        assert missing["items"] == []
        assert missing["pagination"]["total"] == 0
        assert client.get("/v4/aliases?q=x&search_fields=bogus").status_code == 422

    @pytest.mark.component
    def test_resolve_alias_from_real_store(self, tmp_path, fixtures_path):
        """Resolve returns the digest, collection level-2 data, and FHR metadata.

        Uses a real ingested genome (no store mock): initialize_genome
        populates the RefgetStore and apply_fhr writes the sidecar the
        endpoint serves.
        """
        rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        rgc.init(genome_folder=tmp_path / "genomes")
        register_fasta(rgc, fixtures_path)
        digest, _ = rgc.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            description="",
            alias_names=["rCRSd"],
        )
        rgc.genome.apply_fhr(digest, {"genome": "Homo sapiens", "documentation": "A test genome."})

        from refgenie.server.main import create_app

        with TestClient(create_app(refgenie_instance=rgc), raise_server_exceptions=True) as client:
            data = client.get("/v4/aliases/rCRSd").json()
        assert data["alias"] == "rCRSd"
        assert data["digest"] == digest
        assert "names" in data["collection"]
        assert "lengths" in data["collection"]
        assert data["fhr"]["genome"] == "Homo sapiens"

    def test_resolve_alias_not_found(self, client_with_alias):
        """An unknown alias is a 404 (no namespace fallback)."""
        client = client_with_alias
        assert client.get("/v4/aliases/nonexistent").status_code == 404
        assert client.get("/v4/aliases/hg19").status_code == 404

    def test_resolve_alias_collection_not_found(self, client_with_alias, rgc_with_alias):
        """Alias exists but the store has no collection data -> 404 with detail."""
        client, test_rgc = client_with_alias, rgc_with_alias
        mock_router = MagicMock()
        mock_router.get_collection_level2.return_value = None
        test_rgc._store_router = mock_router

        response = client.get("/v4/aliases/rCRSd")
        assert response.status_code == 404
        assert "Collection data not found" in response.json()["detail"]


class TestGenomesEndpoint:
    """GET /v4/genomes in local mode sources aliases from the store, not a SQL join.

    Regression: a select(Genome).join(Alias) inner join drops every genome
    because in local mode aliases live in the RefgetStore.
    """

    def test_list_genomes_populated_with_aliases(self, client_with_alias):
        client = client_with_alias
        data = client.get("/v4/genomes").json()
        genome = next(i for i in data["items"] if i["digest"] == "test_digest_12345")
        assert set(genome["aliases"]) == {"rCRSd", "hg38_mito"}

    def test_list_genomes_filter_by_alias(self, client_with_alias):
        client = client_with_alias
        found = client.get("/v4/genomes?alias=rCRSd").json()
        assert [i["digest"] for i in found["items"]] == ["test_digest_12345"]
        assert client.get("/v4/genomes?alias=nonexistent").json()["items"] == []

    def test_list_genomes_search_by_alias(self, client_with_alias):
        client = client_with_alias
        found = client.get("/v4/genomes?q=mito&search_fields=aliases").json()
        assert "test_digest_12345" in [i["digest"] for i in found["items"]]
        assert client.get("/v4/genomes?q=zzzznomatch&search_fields=aliases").json()["items"] == []

    def test_genome_without_alias_appears(self, client_no_alias):
        """A genome with zero aliases must still be listed, with aliases == []."""
        data = client_no_alias.get("/v4/genomes").json()
        genome = next(i for i in data["items"] if i["digest"] == "no_alias_digest_999")
        assert genome["aliases"] == []


class TestServiceInfo:
    """The three service-info documents: seqcol, root discovery, and DRS."""

    def test_seqcol_service_info_baseline(self, service_info_client):
        """GA4GH identity fields, and refget_store disabled without a store URL."""
        data = service_info_client.get("/seqcol/service-info").json()
        assert data["id"] == "org.refgenie.seqcol"
        assert data["name"] == "Refgenie Sequence Collections"
        assert data["type"]["group"] == "org.ga4gh"
        assert data["type"]["artifact"] == "refget-seqcol"
        assert data["seqcol"]["refget_store"]["enabled"] is False
        assert "url" not in data["seqcol"]["refget_store"]

    def test_root_service_info_identifies_refgenie(self, service_info_client):
        data = service_info_client.get("/service-info").json()
        assert data["id"] == "org.refgenie.api"
        assert data["name"] == "Refgenie"
        assert data["type"]["group"] == "org.refgenie"
        assert data["type"]["artifact"] == "refgenie"
        assert data["version"]

    def test_root_service_info_carries_the_seqcol_bootstrap_shape(self, service_info_client):
        """seqcol.refget_store is what a client reads to learn the store URL.

        Without a store it reports enabled: false, which makes the client refuse
        this server as a genome source loudly instead of registering
        sequence-less genomes.
        """
        seqcol = service_info_client.get("/service-info").json()["seqcol"]
        assert seqcol["url"] == "/seqcol"
        assert seqcol["service_info"] == "/seqcol/service-info"
        assert seqcol["id"] == "org.refgenie.seqcol"
        assert seqcol["refget_store"]["enabled"] is False

    def test_root_service_info_lists_seqcol_mounts(self, service_info_client):
        """The mount is also given as a list, so more can be added without a break."""
        data = service_info_client.get("/service-info").json()
        services = data["seqcol"]["services"]
        assert isinstance(services, list)
        assert len(services) == 1
        assert services[0]["url"] == "/seqcol"
        assert services[0]["refget_store"] == data["seqcol"]["refget_store"]

    def test_root_service_info_is_parsed_by_the_client_bootstrap(self, service_info_client):
        """The client's parser reads this document, returning None when disabled."""
        from refgenie.managers.sources.genomes import store_url_from_service_info

        doc = service_info_client.get("/service-info").json()
        assert store_url_from_service_info(doc) is None

    def test_drs_service_info_still_works(self, service_info_client):
        """Adding a root /service-info must not disturb the DRS one at both mounts."""
        for path in ("/ga4gh/drs/service-info", "/v4/ga4gh/drs/service-info"):
            response = service_info_client.get(path)
            assert response.status_code == 200, path
            assert response.json()["id"] == "org.refgenie.server"
            assert response.json()["type"]["artifact"] == "drs"


GENOME_DIGEST = "test_digest_12345"


def _fasta_spec(*, with_fai=False, stage_modes=()):
    """One seeded fasta asset: the on-disk files plus how it is staged."""
    files = {"test.fa": ">chr1\nACGT\n"}
    if with_fai:
        files["test.fa.fai"] = "chr1\t4\t6\t4\t5\n"
    return {
        "key": "asset",
        "asset_group_name": "fasta",
        "asset_class": "fasta",
        "files": files,
        "stage_modes": stage_modes,
    }


def _server_client(tmp_path, fixtures_path, assets):
    """(client, digests) over a server refgenie seeded with ``assets``."""
    rgc, digests = server_with_asset(
        tmp_path,
        genome_digest=GENOME_DIGEST,
        fixtures_path=fixtures_path,
        assets=assets,
    )
    return make_server_client(rgc), digests


@pytest.fixture
def server_client_minimal(tmp_path, fixtures_path):
    """One genome, no assets."""
    tc, digests = _server_client(tmp_path, fixtures_path, ())
    with tc as client:
        yield client, digests


@pytest.fixture
def server_client_with_asset(tmp_path, fixtures_path):
    """One unstaged fasta asset."""
    tc, digests = _server_client(tmp_path, fixtures_path, [_fasta_spec()])
    with tc as client:
        yield client, digests


@pytest.fixture
def server_client_with_file_asset(tmp_path, fixtures_path):
    """One file-mode-staged fasta asset (with a .fai)."""
    tc, digests = _server_client(
        tmp_path, fixtures_path, [_fasta_spec(with_fai=True, stage_modes=("file",))]
    )
    with tc as client:
        yield client, digests


@pytest.fixture
def server_client_with_archive_asset(tmp_path, fixtures_path):
    """One archive-mode-staged fasta asset."""
    tc, digests = _server_client(
        tmp_path, fixtures_path, [_fasta_spec(stage_modes=("archive",))]
    )
    with tc as client:
        yield client, digests


@pytest.fixture
def core_client_with_asset(tmp_path, fixtures_path):
    """One unstaged fasta asset (with a .fai), for the /v4/assets shape tests."""
    tc, digests = _server_client(tmp_path, fixtures_path, [_fasta_spec(with_fai=True)])
    with tc as client:
        yield client, digests


class TestServerAPIEndpoints:
    """Health and summary endpoints return valid responses."""

    pytestmark = pytest.mark.component

    def test_healthcheck(self, server_client_minimal):
        client, _ = server_client_minimal
        response = client.get("/v4/healthcheck")
        assert response.status_code == 200
        assert response.json()["status"] == "ok"

    def test_summary(self, server_client_minimal):
        client, _ = server_client_minimal
        data = client.get("/v4/summary").json()
        assert "genomes" in data
        assert "asset_groups" in data
        assert "assets" in data

    def test_species_summary(self, server_client_minimal):
        client, _ = server_client_minimal
        assert client.get("/v4/species/summary").status_code == 200


class TestFileServing:
    """/v4/assets/{digest}/files listing and download, with path-traversal + mode guards."""

    pytestmark = pytest.mark.component

    def test_list_asset_files(self, server_client_with_file_asset):
        client, digests = server_client_with_file_asset
        response = client.get(f"/v4/assets/{digests['asset']}/files")
        assert response.status_code == 200, response.text
        data = response.json()
        file_list = data["files"] if isinstance(data, dict) and "files" in data else data
        assert "test.fa" in file_list
        assert "test.fa.fai" in file_list

    def test_list_asset_files_not_found(self, server_client_with_file_asset):
        client, _ = server_client_with_file_asset
        assert client.get("/v4/assets/nonexistent_digest/files").status_code == 404

    def test_download_asset_file(self, server_client_with_file_asset):
        client, digests = server_client_with_file_asset
        response = client.get(f"/v4/assets/{digests['asset']}/files/test.fa")
        assert response.status_code == 200
        assert ">chr1" in response.text

    def test_download_asset_file_path_traversal(self, server_client_with_file_asset):
        """A traversal path must 404, never escape the asset directory."""
        client, digests = server_client_with_file_asset
        response = client.get(f"/v4/assets/{digests['asset']}/files/../../etc/passwd")
        assert response.status_code == 404

    def test_download_asset_file_not_in_listing(self, server_client_with_file_asset):
        client, digests = server_client_with_file_asset
        response = client.get(f"/v4/assets/{digests['asset']}/files/nonexistent.txt")
        assert response.status_code == 404

    def test_download_archive_no_archive_staged(self, server_client_with_file_asset):
        """A file-mode-only asset has no archive to download -> 404."""
        client, digests = server_client_with_file_asset
        assert client.get(f"/v4/archives/{digests['asset']}/download").status_code == 404


class TestArchivesList:
    """GET /v4/archives -- the per-asset archive metadata the UI's tables need."""

    pytestmark = pytest.mark.component

    def test_list_archives(self, server_client_with_archive_asset):
        client, digests = server_client_with_archive_asset
        items = client.get("/v4/archives").json()["items"]
        assert len(items) == 1
        # `digest` must be the asset digest: the UI feeds it into /archives/{digest}/download.
        assert items[0]["digest"] == digests["asset"]
        assert items[0]["asset_digest"] == digests["asset"]
        assert items[0]["size"] == 2048
        assert items[0]["directory_contents"] == ["test.fa"]

    def test_file_mode_not_listed(self, server_client_with_file_asset):
        """File-mode staging is not an archive; the list must exclude it."""
        client, _ = server_client_with_file_asset
        response = client.get("/v4/archives")
        assert response.status_code == 200
        assert response.json()["items"] == []

    @pytest.mark.parametrize(
        "query,expected",
        [
            ("asset_digest={digest}", 1),
            ("asset_digest=" + "0" * 64, 0),
            ("genome_digest=test_digest_12345", 1),
            ("genome_digest=" + "0" * 32, 0),
            ("genome_digest=test_digest_12345&asset_digest={digest}", 1),  # AND: both match
            ("genome_digest=test_digest_12345&asset_digest=" + "0" * 64, 0),  # AND: asset wrong
        ],
    )
    def test_archive_filters(self, server_client_with_archive_asset, query, expected):
        client, digests = server_client_with_archive_asset
        response = client.get(f"/v4/archives?{query.format(digest=digests['asset'])}")
        assert response.status_code == 200
        assert len(response.json()["items"]) == expected


class TestAssetServingModesInResponse:
    """The assets endpoints expose resolved serving_modes so pulls work."""

    pytestmark = pytest.mark.component

    def test_get_asset_includes_resolved_serving_modes(self, core_client_with_asset):
        client, digests = core_client_with_asset
        response = client.get(f"/v4/assets/{digests['asset']}")
        assert response.status_code == 200, response.text
        data = response.json()
        # The fasta asset class ships serving_modes=['file', 'archive']; no override.
        assert data["serving_modes"] == ["file", "archive"]
        assert data["asset_class_name"] == "fasta"
        # Group name and genome digest are embedded so the UI needs no extra fetch.
        assert data["asset_group_name"] == "fasta"
        assert data["genome_digest"] == GENOME_DIGEST

    def test_list_assets_includes_resolved_serving_modes(self, core_client_with_asset):
        """The list path serializes serving_modes without lazy-load errors.

        raise_server_exceptions=True means a lazy-load error during serialization
        surfaces here instead of a clean 200.
        """
        client, _ = core_client_with_asset
        response = client.get("/v4/assets")
        assert response.status_code == 200, response.text
        items = response.json()["items"]
        assert items, "expected at least one asset"
        for item in items:
            assert item["serving_modes"] == ["file", "archive"]
            assert item["asset_class_name"] == "fasta"
            assert item["asset_group_name"] == "fasta"
            assert item["genome_digest"] == GENOME_DIGEST


# --- to_public() / expand=true field completeness ----------------------------


@pytest.mark.component
def test_relationships_expand_returns_the_asset_class_and_serving_modes(tmp_path, fixtures_path):
    """`?expand=true` must return real asset objects, not stripped stubs.

    `to_public` built `asset_class_name`, `serving_modes`, `parents` and
    `children` and then handed them to `AssetPublic`, which declares none of
    them -- pydantic's `extra="ignore"` dropped all four silently.
    """
    test_rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    genome_folder = tmp_path / "genomes"
    test_rgc.init(genome_folder=genome_folder)
    register_fasta(test_rgc, fixtures_path)
    test_rgc.genome.add(digest="d" * 16, description="g", alias_names=["g1"])
    asset_class = test_rgc.asset_class.get("fasta")

    digests = {}
    with Session(test_rgc.database_engine) as session:
        group = AssetGroup(
            name="fasta",
            description="grp",
            genome_digest="d" * 16,
            asset_class_id=asset_class.id,
        )
        session.add(group)
        session.commit()
        session.refresh(group)
        for role in ("parent", "child"):
            adir = genome_folder / "data" / role
            adir.mkdir(parents=True, exist_ok=True)
            (adir / "x.fa").write_text(f">{role}\nACGT\n")
            d = get_dir_digest(adir)
            session.add(
                Asset(
                    digest=d, name=role, description=role, asset_group_id=group.id, path=str(adir)
                )
            )
            session.add(
                AssetName(
                    name=role, asset_group_id=group.id, asset_digest=d, is_default=role == "parent"
                )
            )
            digests[role] = d
        session.commit()
        child = session.get(Asset, digests["child"])
        child.parents.append(session.get(Asset, digests["parent"]))
        session.add(child)
        session.commit()

    with make_server_client(test_rgc) as client:
        body = client.get(f"/v4/relationships/{digests['child']}?expand=true").json()

    assert [p["digest"] for p in body["parents"]] == [digests["parent"]]
    assert body["parents"][0]["asset_class_name"] == "fasta"
    assert body["parents"][0]["serving_modes"]


# --- core endpoints ----------------------------------------------------------
#
# The local-mode app's /v4 router over an empty in-memory DB with a stubbed
# refgenie. These stay `unit`: nothing is built on disk.

_test_engine = make_engine()
SQLModel.metadata.create_all(_test_engine)


@pytest.fixture(scope="module")
def local_client():
    """Test client for the shared /v4 router over the empty in-memory DB.

    This is the REAL local-mode app (``create_app(mode="local")``), not a
    hand-built FastAPI(): the thing that ships is the thing under test. Its
    Refgenie is a stub and its sessions come from the in-memory engine, so
    nothing here touches the caller's real configuration.
    """
    with make_local_client(stub_rgc(), engine=_test_engine) as c:
        yield c


_CORE_ENDPOINTS = [
    "/v4/genomes",
    "/v4/asset_groups",
    "/v4/assets",
    "/v4/asset_classes",
    "/v4/recipes",
    "/v4/configurations",
    "/v4/staged_assets",
]


class TestCorePagination:
    """Every paginated /v4 endpoint shares one envelope and bounds contract."""

    @pytest.mark.parametrize("endpoint", _CORE_ENDPOINTS)
    def test_pagination_envelope(self, local_client, endpoint):
        data = local_client.get(endpoint).json()
        assert isinstance(data["items"], list)
        pagination = data["pagination"]
        assert pagination["offset"] == 0
        assert pagination["limit"] == 100
        assert isinstance(pagination["total"], int)

    @pytest.mark.parametrize("endpoint", _CORE_ENDPOINTS)
    def test_pagination_bounds_and_echo(self, local_client, endpoint):
        assert local_client.get(f"{endpoint}?offset=-1").status_code == 422
        assert local_client.get(f"{endpoint}?limit=0").status_code == 422
        assert local_client.get(f"{endpoint}?limit=1001").status_code == 422
        # max boundary accepted
        assert local_client.get(f"{endpoint}?limit=1000").status_code == 200
        page = local_client.get(f"{endpoint}?offset=5&limit=10").json()
        assert page["pagination"]["offset"] == 5
        assert page["pagination"]["limit"] == 10

    @pytest.mark.parametrize(
        "endpoint,filter_param",
        [
            ("/v4/genomes", "q=nonexistent_digest_12345&search_fields=digest"),
            ("/v4/asset_groups", "q=nonexistent_12345&search_fields=name"),
            ("/v4/assets", "q=nonexistent_asset_12345&search_fields=name"),
            ("/v4/asset_classes", "q=nonexistent_class_12345&search_fields=name"),
            ("/v4/recipes", "q=nonexistent_recipe_12345&search_fields=name"),
            ("/v4/staged_assets", "q=nonexistent_staged_12345&search_fields=asset_digest"),
        ],
    )
    def test_empty_results(self, local_client, endpoint, filter_param):
        data = local_client.get(f"{endpoint}?{filter_param}&offset=0&limit=10").json()
        assert data["items"] == []
        assert data["pagination"]["total"] == 0
        assert data["pagination"]["limit"] == 10


class TestCoreLookups:
    """Single-item lookups 404, accepted search fields 200, bad fields 422."""

    @pytest.mark.parametrize(
        "path",
        [
            "/v4/genomes/nonexistent",
            "/v4/asset_groups/999999",
            "/v4/assets/nonexistent",
            "/v4/assets/nonexistent/files",
            "/v4/asset_classes/999999",
            "/v4/recipes/999999",
            "/v4/configurations/999999",
            "/v4/staged_assets/999",
            "/v4/relationships/nonexistent",
        ],
    )
    def test_get_single_404(self, local_client, path):
        assert local_client.get(path).status_code == 404

    @pytest.mark.parametrize(
        "endpoint,field",
        [
            ("/v4/genomes", "digest"),
            ("/v4/genomes", "aliases"),
            ("/v4/asset_groups", "name"),
            ("/v4/assets", "name"),
            ("/v4/assets", "digest"),
            ("/v4/assets", "path"),
            ("/v4/asset_classes", "name"),
            ("/v4/asset_classes", "version"),
            ("/v4/recipes", "name"),
            ("/v4/recipes", "version"),
            ("/v4/recipes", "description"),
            ("/v4/staged_assets", "asset_digest"),
            ("/v4/staged_assets", "mode"),
        ],
    )
    def test_search_field_accepted(self, local_client, endpoint, field):
        assert (
            local_client.get(f"{endpoint}?q=nonexistent&search_fields={field}").status_code == 200
        )

    def test_invalid_search_fields(self, local_client):
        """Invalid search fields are a 422 that names the offenders and valid set.

        The local app wraps every non-2xx in the ``{"ok": false, "error": ...}``
        envelope (refgenie/server/errors.py), so the message lives at
        ``error.message`` here rather than FastAPI's bare ``detail``.
        """
        r = local_client.get("/v4/genomes?q=test&search_fields=invalid_field")
        assert r.status_code == 422
        assert "Invalid search fields: invalid_field" in r.json()["error"]["message"]
        assert "Valid fields are:" in r.json()["error"]["message"]

        mixed = local_client.get(
            "/v4/genomes?q=test&search_fields=digest,invalid_field,species_name"
        )
        assert "Invalid search fields: invalid_field" in mixed.json()["error"]["message"]

        multi = local_client.get("/v4/genomes?q=test&search_fields=invalid1,invalid2")
        assert "Invalid search fields: invalid1, invalid2" in multi.json()["error"]["message"]


class TestLocalRouteOwnership:
    """Each path is defined by one handler on one router.

    See docs/design-notes.md ("Route ownership across app modes"): the formerly
    duplicated alias and asset-files paths are defined once, on the shared
    router, which both modes include at /v4.
    """

    def test_local_mode_serves_alias_endpoints(self, local_client):
        data = local_client.get("/v4/aliases").json()
        assert "items" in data
        assert "pagination" in data
        # Unknown alias -> 404, not 500 or an empty 200.
        assert local_client.get("/v4/aliases/definitely-not-a-real-alias").status_code == 404


# --- data-channel router -----------------------------------------------------

# A native-layout index like refgenie-registry publishes: flat asset classes,
# nested recipes (<name>/recipe.yaml).
NATIVE_INDEX = IndexFile(
    asset_class=IndexFileSection(dir="asset_classes", files=["bwa_index.yaml"]),
    recipe=IndexFileSection(dir="recipes", files=["bwa_index/recipe.yaml"]),
)


@pytest.fixture
def data_channel_client(monkeypatch, tmp_path, fixtures_path):
    """A TestClient over the real server app, upstream index fetch stubbed."""

    async def fake_fetch(self):
        self.index_data = NATIVE_INDEX

    # The default config ships a 'registry' channel pointing at the Pages URL;
    # stub the network fetch so tests are hermetic but redirects still use the
    # real configured index_url.
    monkeypatch.setattr(DataChannel, "fetch_index_yaml", fake_fetch)

    with make_server_client(make_server_rgc(tmp_path, fixtures_path)) as client:
        yield client


def test_lists_configured_channels(data_channel_client):
    r = data_channel_client.get("/v4/data_channel/")
    assert r.status_code == 200
    assert "registry" in r.json()


def test_channel_index_served_verbatim_native(data_channel_client):
    r = data_channel_client.get("/v4/data_channel/registry/index.yaml")
    assert r.status_code == 200
    idx = yaml.safe_load(r.text)
    # Native structure preserved: dirs + nested recipe path, no rewriting.
    assert idx["asset_class"]["dir"] == "asset_classes"
    assert idx["recipe"]["dir"] == "recipes"
    assert idx["recipe"]["files"] == ["bwa_index/recipe.yaml"]
    assert idx["asset_class"]["files"] == ["bwa_index.yaml"]


def test_nested_recipe_path_redirects_to_upstream(data_channel_client):
    r = data_channel_client.get(
        "/v4/data_channel/registry/recipes/bwa_index/recipe.yaml",
        follow_redirects=False,
    )
    assert r.status_code == 302
    # urljoin against the configured index_url replaces the trailing index.yaml.
    assert r.headers["location"].endswith("/recipes/bwa_index/recipe.yaml")
    assert "__" not in r.headers["location"]


def test_asset_class_path_redirects_to_upstream(data_channel_client):
    r = data_channel_client.get(
        "/v4/data_channel/registry/asset_classes/bwa_index.yaml",
        follow_redirects=False,
    )
    assert r.status_code == 302
    assert r.headers["location"].endswith("/asset_classes/bwa_index.yaml")


def test_unknown_channel_is_404(data_channel_client):
    assert data_channel_client.get("/v4/data_channel/nope/index.yaml").status_code == 404
    assert (
        data_channel_client.get("/v4/data_channel/nope/recipes/x/recipe.yaml").status_code == 404
    )


# --- Route-registration guards ----------------------------------------------
#
# The server app mounts several routers at the same prefixes ("" and "/v4").
# FastAPI resolves collisions first-match-wins, so a duplicate route silently
# becomes dead code -- but the published OpenAPI still advertises both, and
# clients generated from it are wrong. These guard against that.


def _normalize(path):
    """Drop Starlette path converters ({p:path} -> {p}) as OpenAPI does."""
    return re.sub(r"\{([^}:]+):[^}]+\}", r"{\1}", path)


def _effective_api_routes(app):
    """Yield an effective view of every APIRoute reachable from `app`.

    As of fastapi 0.137 `include_router` no longer copies child routes onto
    `app.routes`; it appends one opaque `_IncludedRouter` and computes the
    effective (prefixed) routes lazily. `iter_route_contexts` is the walker
    `get_openapi` itself uses; each yielded `RouteContext` exposes the effective
    path, methods and operation id, delegating to the underlying route for
    routes declared directly on the app (which need no prefixing).

    The duplicate-registration bug these guards exist to catch is invisible in
    the public OpenAPI schema -- including a router twice yields a byte-identical
    spec -- so there is no public alternative to this walker. If the running
    FastAPI does not expose it, skip rather than assert on a partial view.
    """
    if iter_route_contexts is None:
        pytest.skip(
            "FastAPI does not expose routing.iter_route_contexts; cannot walk "
            "effective routes to guard against duplicate registration"
        )
    for ctx in iter_route_contexts(app.router.routes):
        if isinstance(ctx.original_route, APIRoute):
            yield ctx


def _route_keys(app):
    """Yield (method, path) for every APIRoute on `app`."""
    for route in _effective_api_routes(app):
        for method in sorted(route.methods or []):
            yield method, _normalize(route.path)


def _assert_no_duplicate_routes(app, label):
    counts = Counter(_route_keys(app))
    duplicates = {key: n for key, n in counts.items() if n > 1}
    assert not duplicates, f"{label}: route(s) registered more than once: " + ", ".join(
        f"{m} {p} (x{n})" for (m, p), n in sorted(duplicates.items())
    )


def _assert_unique_operation_ids(app, label):
    """Duplicate operationIds break generated clients even when paths differ."""
    seen = defaultdict(list)
    for route in _effective_api_routes(app):
        if route.include_in_schema:
            seen[route.operation_id or route.unique_id].append(route.path)
    duplicates = {op: paths for op, paths in seen.items() if len(paths) > 1}
    assert not duplicates, f"{label}: duplicate operationIds: {duplicates}"


class TestRouteUniqueness:
    """No HTTP route or operationId may be registered twice."""

    def test_server_app_has_no_duplicate_routes(self, server_app):
        _assert_no_duplicate_routes(server_app, "server app")

    def test_server_app_has_unique_operation_ids(self, server_app):
        _assert_unique_operation_ids(server_app, "server app")

    def test_server_openapi_paths_are_unambiguous(self, server_app):
        """Every documented path+method must map to exactly one registered route."""
        spec = server_app.openapi()
        registered = Counter(_route_keys(server_app))
        for path, operations in spec["paths"].items():
            for method in operations:
                key = (method.upper(), path)
                assert registered[key] == 1, (
                    f"OpenAPI advertises {method.upper()} {path} but "
                    f"{registered[key]} route(s) are registered for it"
                )


class TestEndpointHitCounter:
    """The hit counter must see routes added via include_router.

    It used to hand-roll route matching over `app.router.routes` before dispatch.
    fastapi 0.137 made that silently record nothing for every included route --
    breaking archive download counts -- without failing any test.
    """

    @staticmethod
    def _app_with_counter():
        from fastapi import APIRouter, FastAPI

        from refgenie.server.stats import EndpointCollector, EndpointHitCounterMiddleware

        inner = APIRouter()

        @inner.get("/archive/{asset_digest}", name="download_archive")
        def download_archive(asset_digest: str):
            return {"asset": asset_digest}

        outer = APIRouter()
        outer.include_router(inner)

        app = FastAPI()
        app.include_router(outer, prefix="/v4")
        app.state.endpoint_hits_collector = EndpointCollector()
        return app, EndpointHitCounterMiddleware(app)

    def test_records_hit_for_included_route_with_path_param(self):
        app, wrapped = self._app_with_counter()
        with TestClient(wrapped) as client:
            assert client.get("/v4/archive/abc123").status_code == 200

        stats = app.state.endpoint_hits_collector.get_stats("download_archive")
        assert stats, "no hit recorded for a prefixed include_router route"
        assert list(stats) == [frozenset({("asset_digest", "abc123")})]
        assert app.state.endpoint_hits_collector.get_total_hits("download_archive") == 1

    def test_records_nothing_for_unmatched_path(self):
        app, wrapped = self._app_with_counter()
        with TestClient(wrapped) as client:
            assert client.get("/v4/nonexistent").status_code == 404

        assert app.state.endpoint_hits_collector.get_stats("download_archive") == {}


def test_get_db_session_honours_the_refgenie_override(tmp_path):
    """A handler taking both `rgc` and `session` must get one database.

    `create_app(refgenie_instance=X)` installs
    ``app.dependency_overrides[get_refgenie] = lambda: X``. That override is
    only reachable if `get_db_session` declares `get_refgenie` as a dependency
    rather than calling it.
    """
    from fastapi import Depends, FastAPI
    from fastapi.testclient import TestClient
    from sqlmodel import Session

    from refgenie.server import dependencies as dash_dependencies

    def _rgc():
        return Refgenie(database_engine=make_engine(), suppress_migrations=True)

    singleton, override = _rgc(), _rgc()
    assert singleton.database_engine is not override.database_engine

    app = FastAPI()

    @app.get("/engine")
    def read_engine(session: Session = Depends(dash_dependencies.get_db_session)):
        return {"engine": str(id(session.get_bind()))}

    app.dependency_overrides[dash_dependencies.get_refgenie] = lambda: override

    saved = dash_dependencies._refgenie_instance
    dash_dependencies._refgenie_instance = singleton
    try:
        with TestClient(app) as client:
            body = client.get("/engine").json()
    finally:
        dash_dependencies._refgenie_instance = saved

    assert body["engine"] == str(id(override.database_engine))
