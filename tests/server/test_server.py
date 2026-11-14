"""The server app over HTTP: v4 endpoints, service-info, file serving, the core
router and the data-channel mirror.
"""

from collections import Counter
from unittest.mock import MagicMock
from refgenie.exceptions import MissingGenomeError

import pytest
import yaml

from tests.support.factories import (
    fake_digest,
    make_engine,
    make_server_rgc,
    register_fasta,
    seed_search_genomes,
    seed_species,
)
from tests.support.server import (
    make_local_client,
    make_server_app,
    make_server_client,
    requires_dash,
    requires_server,
    serve_refgenie,
    server_with_asset,
)

TEST_DIGEST = fake_digest("test_digest_12345")
NO_ALIAS_DIGEST = fake_digest("no_alias_digest_999")

# Server extras are optional; skip the whole module if any is missing.
requires_server()
requires_dash()

from fastapi.testclient import TestClient  # noqa: E402  (must follow the extras guard)
from pydantic import ValidationError  # noqa: E402
from sqlmodel import Session  # noqa: E402

from refgenie import Refgenie  # noqa: E402
from refgenie.db.tables import (  # noqa: E402
    Asset,
    AssetGroup,
    AssetName,
    Genome,
    SeekKey,
    SeekKeyType,
    StoreType,
)
from refgenie.managers.sources.manager import IndexFile, IndexFileSection  # noqa: E402
from refgenie.server.routers.data_channel import resolve_channel_file_url  # noqa: E402
from refgenie.server.schemas import DataChannel  # noqa: E402
from refgenie.utils.digest import get_dir_digest  # noqa: E402

# --- v4 server endpoints -----------------------------------------------------


@pytest.fixture
def rgc_with_alias(tmp_path, fixtures_path):
    """A catalog holding one genome (TEST_DIGEST) with two aliases."""
    return make_server_rgc(tmp_path, fixtures_path, genomes=[(TEST_DIGEST, ["rCRSd", "hg38_mito"])])


@pytest.fixture
def client_with_alias(rgc_with_alias):
    """Client over ``rgc_with_alias``."""
    with make_server_client(rgc_with_alias) as client:
        yield client


@pytest.fixture
def client_no_alias(tmp_path, fixtures_path):
    """Client over a genome with no aliases."""
    rgc = make_server_rgc(tmp_path, fixtures_path, genomes=[(NO_ALIAS_DIGEST, [])])
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
        names = [item["name"] for item in data["results"]]
        assert "rCRSd" in names
        assert "hg38_mito" in names

    def test_list_aliases_filter_by_genome(self, client_with_alias):
        client = client_with_alias
        assert len(client.get(f"/v4/aliases?genome_digest={TEST_DIGEST}").json()["results"]) == 2
        assert len(client.get("/v4/aliases?genome_digest=nonexistent").json()["results"]) == 0

    def test_list_aliases_is_paginated(self, client_with_alias):
        """The client's page loop needs `pagination` and honored page/page_size.

        An unpaginated list makes ServerSource.get_paginated repeat forever once a
        server has more aliases than one page.
        """
        client = client_with_alias
        assert client.get("/v4/aliases").json()["pagination"] == {
            "page": 0,
            "page_size": 100,
            "total": 2,
        }
        page = client.get("/v4/aliases?page=1&page_size=1").json()
        assert page["pagination"] == {"page": 1, "page_size": 1, "total": 2}
        assert [i["name"] for i in page["results"]] == ["rCRSd"]

    @pytest.mark.parametrize(
        ["path", "status"],
        [
            ("/v4/aliases?page=-1", 422),
            ("/v4/aliases?page_size=0", 422),
            ("/v4/aliases?page_size=1001", 400),
            ("/v4/aliases?page=2&page_size=1", 400),
        ],
    )
    def test_list_aliases_rejects_bad_pagination(self, client_with_alias, path, status):
        client = client_with_alias
        assert client.get(path).status_code == status

    def test_list_aliases_filter_by_name_and_search(self, client_with_alias):
        client = client_with_alias
        assert [i["name"] for i in client.get("/v4/aliases?name=rCRSd").json()["results"]] == [
            "rCRSd"
        ]
        found = client.get("/v4/aliases?q=mito&search_fields=name").json()
        assert [i["name"] for i in found["results"]] == ["hg38_mito"]
        missing = client.get("/v4/aliases?q=zzzznomatch").json()
        assert missing["results"] == []
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
        rgc.database.init(genome_folder=tmp_path / "genomes")
        register_fasta(rgc, fixtures_path)
        digest, _ = rgc.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            description="",
            alias_names=["rCRSd"],
        )
        rgc.genome.apply_fhr(
            digest,
            {
                "genome": "rCRSd",
                "version": "1",
                "taxon": {"name": "Homo sapiens", "uri": "https://identifiers.org/taxonomy:9606"},
                "documentation": "A test genome.",
            },
        )

        from refgenie.server.main import create_app

        with TestClient(create_app(refgenie_instance=rgc), raise_server_exceptions=True) as client:
            data = client.get("/v4/aliases/rCRSd").json()
        assert data["alias"] == "rCRSd"
        assert data["digest"] == digest
        assert "names" in data["collection"]
        assert "lengths" in data["collection"]
        assert data["fhr"]["genome"] == "rCRSd"

    def test_resolve_alias_not_found(self, client_with_alias):
        """An unknown alias is a 404 (no namespace fallback)."""
        client = client_with_alias
        assert client.get("/v4/aliases/nonexistent").status_code == 404
        assert client.get("/v4/aliases/hg19").status_code == 404

    def test_resolve_alias_collection_not_found(self, client_with_alias, rgc_with_alias):
        """Alias exists but the store has no collection data -> 404 with detail."""
        client, test_rgc = client_with_alias, rgc_with_alias
        mock_router = MagicMock()
        mock_router.serves.return_value = True
        mock_router.get_collection_level2.side_effect = MissingGenomeError(TEST_DIGEST)
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
        genome = next(i for i in data["results"] if i["digest"] == TEST_DIGEST)
        assert set(genome["aliases"]) == {"rCRSd", "hg38_mito"}

    def test_list_genomes_filter_by_alias(self, client_with_alias):
        client = client_with_alias
        found = client.get("/v4/genomes?alias=rCRSd").json()
        assert [i["digest"] for i in found["results"]] == [TEST_DIGEST]
        assert client.get("/v4/genomes?alias=nonexistent").json()["results"] == []

    def test_list_genomes_search_by_alias(self, client_with_alias):
        client = client_with_alias
        found = client.get("/v4/genomes?q=mito&search_fields=aliases").json()
        assert TEST_DIGEST in [i["digest"] for i in found["results"]]
        assert client.get("/v4/genomes?q=zzzznomatch&search_fields=aliases").json()["results"] == []

    def test_get_genome_releases_db_connection_before_store_lookup(
        self, client_with_alias, rgc_with_alias
    ):
        """GET /v4/genomes/{digest} closes its session before the FHR lookup.

        The first lookup of a genome in a remote store is a network fetch; a
        database connection held idle across it can be dropped underneath.
        """
        from sqlalchemy import event

        engine = rgc_with_alias.database_engine
        outstanding = Counter()
        event.listen(engine, "checkout", lambda *a: outstanding.update(["open"]))
        event.listen(engine, "checkin", lambda *a: outstanding.subtract(["open"]))
        seen = []

        class _Router:
            def get_fhr_metadata(self, digest, store_name=None):
                seen.append(outstanding["open"])
                return None

            def serves(self, store_name):
                return True

            def store_configs(self):
                return []

        # Read the store's aliases now; the stub router below has no store.
        rgc_with_alias.alias.get_for_genome(TEST_DIGEST)
        rgc_with_alias._store_router = _Router()
        response = client_with_alias.get(f"/v4/genomes/{TEST_DIGEST}")
        assert response.status_code == 200
        assert response.json()["digest"] == TEST_DIGEST
        assert seen == [0]

    @pytest.mark.parametrize(
        "store_name, store_type, url, expected",
        [
            (
                "jungle",
                StoreType.remote,
                "https://s3.test/refget-store/jungle",
                "https://s3.test/refget-store/jungle/",
            ),
            ("jungle", StoreType.remote, "https://s3.test/jungle/", "https://s3.test/jungle/"),
            ("disk", StoreType.on_disk, "/data/disk.refget_store", None),
            (None, None, None, None),
        ],
    )
    def test_get_genome_reports_remote_store_url(
        self, client_with_alias, rgc_with_alias, store_name, store_type, url, expected
    ):
        """refget_store_url is the remote store's base URL, else null."""
        digest = fake_digest("stored_genome")
        if store_name:
            rgc_with_alias.store.add(store_name, url, store_type, priority=10)
        rgc_with_alias.genome.add(digest, "stored genome", [], store_name=store_name)
        rows = rgc_with_alias.store.list_all()

        class _Router:
            def get_fhr_metadata(self, digest, store_name=None):
                return None

            def serves(self, store_name):
                return True

            def store_configs(self):
                return rows

        rgc_with_alias.alias.get_for_genome(digest)
        rgc_with_alias._store_router = _Router()
        response = client_with_alias.get(f"/v4/genomes/{digest}")
        assert response.status_code == 200
        assert response.json()["refget_store_url"] == expected

    def test_genome_in_an_unserved_store_is_hidden(self, client_with_alias, rgc_with_alias):
        """A genome whose store_name this server does not serve: 404, and not listed."""
        digest = fake_digest("plant_genome")
        rgc_with_alias.store.add("jungle", "https://s3.test/jungle/", StoreType.remote, priority=10)
        rgc_with_alias.store.add(
            "plantref", "https://s3.test/plant/", StoreType.remote, priority=50
        )
        rgc_with_alias.genome.add(digest, "a plant", [], store_name="plantref")

        class _Router:
            def serves(self, store_name):
                return store_name in (None, "jungle")

            def served_store_names(self):
                return ["jungle"], True

        # Read the store's aliases now; the stub router below has no store.
        rgc_with_alias.alias.get_for_genome(digest)
        rgc_with_alias._store_router = _Router()
        response = client_with_alias.get(f"/v4/genomes/{digest}")
        assert response.status_code == 404
        assert (
            "'plantref', which this server does not serve" in response.json()["detail"]["message"]
        )
        assert client_with_alias.get(f"/v4/genomes/{digest}/sequences").status_code == 404
        listed = [g["digest"] for g in client_with_alias.get("/v4/genomes").json()["results"]]
        assert digest not in listed
        assert TEST_DIGEST in listed
        found = client_with_alias.get("/v4/genomes?q=plant").json()
        assert digest not in [g["digest"] for g in found["results"]]

    def test_genome_without_alias_appears(self, client_no_alias):
        """A genome with zero aliases must still be listed, with aliases == []."""
        data = client_no_alias.get("/v4/genomes").json()
        genome = next(i for i in data["results"] if i["digest"] == NO_ALIAS_DIGEST)
        assert genome["aliases"] == []


class TestGenomePaging:
    """GET /v4/genomes pages in a fixed order: most asset groups first."""

    DIGESTS = [f"{i:032d}" for i in range(7)]
    SPECIES = ["Mus musculus", None, "Homo sapiens", "Homo sapiens", None, "Danio rerio", None]

    @pytest.fixture
    def client(self, refgenie_minimal):
        rgc = refgenie_minimal
        register_fasta(rgc)
        for digest, species in zip(self.DIGESTS, self.SPECIES):
            rgc.genome.add(digest=digest, description="g", alias_names=[], species_name=species)
        fasta_class_id = rgc.asset_class.get("fasta").id
        with Session(rgc.database_engine) as session:
            # Genome 4 (no species) has two asset groups and genome 0 has one, so
            # the count, not the species, puts 4 first.
            for index, name in ((4, "fasta"), (4, "fasta_copy"), (0, "fasta")):
                session.add(
                    AssetGroup(
                        name=name,
                        description="g",
                        genome_digest=self.DIGESTS[index],
                        asset_class_id=fasta_class_id,
                    )
                )
            session.commit()
        rgc.genome.apply_fhr_columns(self.DIGESTS[2], {"taxon": {"uri": "taxonomy:9606"}})
        rgc.genome.apply_fhr_columns(self.DIGESTS[1], {"genomeSynonym": ["maize", "corn"]})
        with make_local_client(rgc) as client:
            yield client

    def test_walking_every_page_returns_each_genome_once(self, client):
        seen, page = [], 0
        while True:
            body = client.get(f"/v4/genomes?page={page}&page_size=2").json()
            assert body["pagination"] == {"page": page, "page_size": 2, "total": 7}
            seen += [g["digest"] for g in body["results"]]
            if (page + 1) * 2 >= body["pagination"]["total"]:
                break
            page += 1
        assert sorted(seen) == sorted(self.DIGESTS)
        assert len(seen) == len(set(seen))
        # Most asset groups first, then by species (unnamed last), then by digest.
        assert seen[:4] == [self.DIGESTS[i] for i in (4, 0, 5, 2)]
        response = client.get(f"/v4/genomes?page={page + 1}&page_size=2")
        assert response.status_code == 400
        assert f"past the last page ({page})" in response.json()["error"]["message"]

    def test_has_assets_and_taxon_id_filter(self, client):
        def digests(query):
            return [g["digest"] for g in client.get(f"/v4/genomes?{query}").json()["results"]]

        assert digests("has_assets=true") == [self.DIGESTS[4], self.DIGESTS[0]]
        assert self.DIGESTS[4] not in digests("has_assets=false")
        assert len(digests("has_assets=false")) == 5
        assert digests("taxon_id=9606") == [self.DIGESTS[2]]
        assert digests("taxon_id=9606&has_assets=true") == []

    def test_synonyms_are_searched_and_returned(self, client):
        def results(query):
            return client.get(f"/v4/genomes?{query}").json()["results"]

        for query in (
            "q=maize",
            "q=maize&operator=eq",
            "q=mai&operator=starts_with",
            "q=orn&operator=ends_with",
            "q=corn&search_fields=synonyms",
        ):
            assert [g["digest"] for g in results(query)] == [self.DIGESTS[1]], query
        assert results("q=mai&operator=eq") == []
        assert results("q=maize")[0]["synonyms"] == ["maize", "corn"]
        detail = client.get(f"/v4/genomes/{self.DIGESTS[1]}").json()
        assert detail["synonyms"] == ["maize", "corn"]


class TestGenomeSpeciesNames:
    """GET /v4/genomes searches the species list and returns its names."""

    HUMAN, MOUSE, BARE_MOUSE = "a" * 32, "b" * 32, "c" * 32
    SPECIES = {
        "taxon_id": 10090,
        "scientific_name": "Mus musculus",
        "common_name": "house mouse",
        "synonyms": ["Hausmaus"],
    }

    @pytest.fixture
    def client(self, refgenie_minimal):
        rgc = seed_search_genomes(refgenie_minimal)
        # A second mouse genome, with no species name of its own.
        with Session(rgc.database_engine) as session:
            session.add(Genome(digest=self.BARE_MOUSE, description="", taxon_id=10090))
            session.commit()
        seed_species(rgc, 10090, "Mus musculus", "house mouse", ["Hausmaus"])
        with make_local_client(rgc) as client:
            yield client

    def _rows(self, client, query):
        return {g["digest"]: g for g in client.get(f"/v4/genomes?{query}").json()["results"]}

    def test_species_list_names_are_searched(self, client):
        for query in ("q=Hausmaus", "q=Hausmaus&search_fields=species"):
            assert sorted(self._rows(client, query)) == [self.MOUSE, self.BARE_MOUSE], query

    def test_common_name_and_species_come_from_the_species_list(self, client):
        rows = self._rows(client, "")
        for digest in (self.MOUSE, self.BARE_MOUSE):
            assert rows[digest]["common_name"] == "house mouse"
            assert rows[digest]["species"] == self.SPECIES
        detail = client.get(f"/v4/genomes/{self.BARE_MOUSE}").json()
        assert detail["common_name"] == "house mouse"
        assert detail["species"] == self.SPECIES

    def test_a_genome_with_no_species_row_has_no_common_name(self, client):
        human = self._rows(client, "")[self.HUMAN]
        assert human["species"] is None
        assert human["common_name"] is None
        assert client.get(f"/v4/genomes/{self.HUMAN}").json()["common_name"] is None


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

    def test_service_info_documents_point_at_the_current_docs(self, service_info_client):
        root = service_info_client.get("/service-info").json()
        drs = service_info_client.get("/ga4gh/drs/v1/service-info").json()
        docs = "https://docs.refgenie.org/refgenie/"
        assert root["documentationUrl"] == drs["documentationUrl"] == docs
        assert root["refgenie"]["links"]["docs"] == docs
        assert root["organization"] == drs["organization"]

    def test_drs_service_info_still_works(self, service_info_client):
        """Adding a root /service-info must not disturb the DRS one."""
        response = service_info_client.get("/ga4gh/drs/v1/service-info")
        assert response.status_code == 200
        assert response.json()["id"] == "org.refgenie.server"
        assert response.json()["type"]["artifact"] == "drs"


GENOME_DIGEST = TEST_DIGEST


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
    tc, digests = _server_client(tmp_path, fixtures_path, [_fasta_spec(stage_modes=("archive",))])
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
        items = client.get("/v4/archives").json()["results"]
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
        assert response.json()["results"] == []

    @pytest.mark.parametrize(
        "query,expected",
        [
            ("asset_digest={digest}", 1),
            ("asset_digest=" + "0" * 64, 0),
            (f"genome_digest={TEST_DIGEST}", 1),
            ("genome_digest=" + "0" * 32, 0),
            (f"genome_digest={TEST_DIGEST}&asset_digest={{digest}}", 1),  # AND: both match
            (f"genome_digest={TEST_DIGEST}&asset_digest=" + "0" * 64, 0),  # AND: asset wrong
        ],
    )
    def test_archive_filters(self, server_client_with_archive_asset, query, expected):
        client, digests = server_client_with_archive_asset
        response = client.get(f"/v4/archives?{query.format(digest=digests['asset'])}")
        assert response.status_code == 200
        assert len(response.json()["results"]) == expected


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
        items = response.json()["results"]
        assert items, "expected at least one asset"
        for item in items:
            assert item["serving_modes"] == ["file", "archive"]
            assert item["asset_class_name"] == "fasta"
            assert item["asset_group_name"] == "fasta"
            assert item["genome_digest"] == GENOME_DIGEST


class TestListAssetsFilters:
    """/v4/assets filters combine, and bad input is a 4xx, never a 500."""

    pytestmark = pytest.mark.component

    @pytest.mark.parametrize(
        "query",
        [
            "asset_group_name=fasta&genome_digest={genome}",
            "asset_group_id={group_id}&genome_digest={genome}",
            "asset_group_id={group_id}&asset_group_name=fasta",
            "asset_group_id={group_id}&asset_group_name=fasta&genome_digest={genome}",
            "name={name}&asset_group_name=fasta&genome_digest={genome}",
        ],
    )
    def test_combined_filters_find_the_asset(self, core_client_with_asset, query):
        client, digests = core_client_with_asset
        asset = client.get(f"/v4/assets/{digests['asset']}").json()
        query = query.format(
            genome=GENOME_DIGEST, group_id=asset["asset_group_id"], name=asset["name"]
        )
        response = client.get(f"/v4/assets?{query}")
        assert response.status_code == 200, response.text
        assert [item["digest"] for item in response.json()["results"]] == [digests["asset"]]

    def test_recipe_filter_combines_with_name_and_genome(self, core_client_with_asset):
        """A bare join(Recipe) after assetname was ambiguous (two FK paths)."""
        client, _ = core_client_with_asset
        response = client.get(
            f"/v4/assets?name=fasta&recipe_name=fasta&genome_digest={GENOME_DIGEST}"
        )
        assert response.status_code == 200, response.text

    @pytest.mark.parametrize(
        "query",
        [
            "asset_group_id=2147483648",
            "asset_group_id=99999999999999999999",
            "asset_group_id=-1",
            "page=99999999999999999999",
            "q=%00",
            "name=%00",
            "asset_group_name=a%00",
            "genome_digest=%00",
            "recipe_name=%00",
        ],
    )
    def test_invalid_input_is_422(self, core_client_with_asset, query):
        client, _ = core_client_with_asset
        assert client.get(f"/v4/assets?{query}").status_code == 422


_ID_ROUTES = [
    "/v4/asset_groups/{id}",
    "/v4/asset_classes/{id}",
    "/v4/recipes/{id}",
    "/v4/configurations/{id}",
    "/v4/staged_assets/{id}",
    "/v4/asset_groups?asset_group_id={id}",
    "/v4/assets?asset_group_id={id}",
]

_NUL_REQUESTS = [
    "/v4/genomes?digest=%00",
    "/v4/genomes?alias=%00",
    "/v4/genomes?q=%00",
    "/v4/asset_groups?asset_class=%00",
    "/v4/asset_classes?name=%00",
    "/v4/asset_classes?version=%00",
    "/v4/recipes?output_asset_class=%00",
    "/v4/staged_assets?asset_digest=%00",
    "/v4/staged_assets?mode=%00",
    "/v4/aliases?name=%00",
    "/v4/archives?digest=%00",
    "/v4/archives?genome_digest=%00",
    "/v4/assets?name=a%00b",
    "/v4/genomes/a%00",
    "/v4/aliases/a%00",
    "/v4/assets/a%00",
    "/v4/relationships/a%00",
    "/v4/archives/a%00/download",
    "/v4/ga4gh/drs/objects/a%00",
]


class TestBadParamsAreClientErrors:
    """Values the database cannot hold are a 422 on every route, never a 500.

    On api.refgenie.org an id past PostgreSQL's INTEGER range (every
    ``/v4/<table>/{id}`` route, and ``asset_group_id``) and a NUL byte in a text
    parameter (nearly every list filter) reached the database and came back as
    a 500. Ids are now typed ``DbId`` and the app rejects NUL bytes in any path
    or query parameter before a handler runs.
    """

    pytestmark = pytest.mark.component

    @pytest.mark.parametrize("route", _ID_ROUTES)
    @pytest.mark.parametrize("value", ["2147483648", "-1", "99999999999999999999"])
    def test_out_of_range_id_is_422(self, server_client_minimal, route, value):
        client, _ = server_client_minimal
        assert client.get(route.format(id=value)).status_code == 422

    @pytest.mark.parametrize("route", _ID_ROUTES)
    def test_largest_id_is_accepted(self, server_client_minimal, route):
        client, _ = server_client_minimal
        assert client.get(route.format(id="2147483647")).status_code in (200, 404)

    @pytest.mark.parametrize("request_path", _NUL_REQUESTS)
    def test_nul_byte_is_422(self, server_client_minimal, request_path):
        client, _ = server_client_minimal
        response = client.get(request_path)
        assert response.status_code == 422, response.text
        assert "NUL" in response.text


def test_local_app_rejects_nul_byte_in_its_error_envelope(local_client):
    """The local app gets the same check, wrapped in its error envelope."""
    response = local_client.get("/v4/genomes?q=%00")
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "validation_error"


# --- to_public() / expand=true field completeness ----------------------------


def _rgc_with_linked_assets(tmp_path, fixtures_path):
    """A server Refgenie with a ``parent`` -> ``child`` asset link.

    Both assets carry a seek key, so an expanded relationship has nested ORM
    rows to serialize. Returns ``(rgc, {"parent": digest, "child": digest})``.
    """
    test_rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    genome_folder = tmp_path / "genomes"
    test_rgc.database.init(genome_folder=genome_folder)
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
            asset = Asset(
                digest=d, name=role, description=role, asset_group_id=group.id, path=str(adir)
            )
            session.add(asset)
            session.add(
                AssetName(
                    name=role, asset_group_id=group.id, asset_digest=d, is_default=role == "parent"
                )
            )
            # Size is computed from the file on insert.
            session.add(
                SeekKey(
                    name="fasta",
                    value="x.fa",
                    description=f"{role} fasta",
                    type=SeekKeyType.file,
                    asset=asset,
                )
            )
            digests[role] = d
        session.commit()
        child = session.get(Asset, digests["child"])
        child.parents.append(session.get(Asset, digests["parent"]))
        session.add(child)
        session.commit()
    return test_rgc, digests


@pytest.mark.component
def test_relationships_expand_returns_the_asset_class_and_serving_modes(tmp_path, fixtures_path):
    """`?expand=true` must return real asset objects, not stripped stubs.

    `to_public` built `asset_class_name`, `serving_modes`, `parents` and
    `children` and then handed them to `AssetPublic`, which declares none of
    them -- pydantic's `extra="ignore"` dropped all four silently.
    """
    test_rgc, digests = _rgc_with_linked_assets(tmp_path, fixtures_path)

    with make_server_client(test_rgc) as client:
        body = client.get(f"/v4/relationships/{digests['child']}?expand=true").json()

    assert [p["digest"] for p in body["parents"]] == [digests["parent"]]
    assert body["parents"][0]["asset_class_name"] == "fasta"
    assert body["parents"][0]["serving_modes"]


@pytest.mark.component
def test_relationships_expand_serializes_nested_seek_keys(tmp_path, fixtures_path):
    """`?expand=true` on an asset whose relatives have seek keys returns 200.

    The handler validated each relative with a bare ``model_validate``, so the
    nested ``SeekKey`` ORM rows reached ``SeekKeyResponse`` without
    ``from_attributes`` and the endpoint raised a 500 for any asset with
    parents. ``refgenie pull`` calls this endpoint with ``expand=true``.
    """
    test_rgc, digests = _rgc_with_linked_assets(tmp_path, fixtures_path)
    expected_seek_keys = lambda role: [  # noqa: E731
        {
            "name": "fasta",
            "value": "x.fa",
            "description": f"{role} fasta",
            "type": "file",
            "size": len(f">{role}\nACGT\n"),
        }
    ]

    with make_server_client(test_rgc) as client:
        up = client.get(f"/v4/relationships/{digests['child']}?expand=true")
        down = client.get(f"/v4/relationships/{digests['parent']}?expand=true")

    assert up.status_code == 200, up.text
    assert down.status_code == 200, down.text
    up, down = up.json(), down.json()

    assert up["asset_digest"] == digests["child"]
    assert up["children"] == []
    (parent,) = up["parents"]
    assert parent["digest"] == digests["parent"]
    assert parent["seek_keys"] == expected_seek_keys("parent")
    assert [n["name"] for n in parent["names"]] == ["parent"]
    assert parent["asset_group_name"] == "fasta"
    assert parent["genome_digest"] == "d" * 16

    assert down["parents"] == []
    (child,) = down["children"]
    assert child["digest"] == digests["child"]
    assert child["seek_keys"] == expected_seek_keys("child")
    assert [n["name"] for n in child["names"]] == ["child"]


# --- core endpoints ----------------------------------------------------------
#
# The local-mode app's /v4 router over an empty in-memory DB. These stay
# `unit`: nothing is built on disk.


@pytest.fixture(scope="module")
def local_client(tmp_path_factory):
    """Test client for the shared /v4 router over an empty in-memory DB.

    This is the REAL local-mode app (``create_app(mode="local")``), not a
    hand-built FastAPI(): the thing that ships is the thing under test. Its
    Refgenie is initialized into a temporary folder, so nothing here touches
    the caller's real configuration.
    """
    root = tmp_path_factory.mktemp("local_client")
    rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rgc.database.init(genome_folder=root / "genomes", genome_stage_folder=root / "archives")
    with make_local_client(rgc) as c:
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
        assert isinstance(data["results"], list)
        pagination = data["pagination"]
        assert pagination["page"] == 0
        assert pagination["page_size"] == 100
        assert isinstance(pagination["total"], int)

    @pytest.mark.parametrize("endpoint", _CORE_ENDPOINTS)
    def test_pagination_bounds_and_echo(self, local_client, endpoint):
        assert local_client.get(f"{endpoint}?page=-1").status_code == 422
        assert local_client.get(f"{endpoint}?page_size=0").status_code == 422
        # Over the cap, and past the end, are both 400s that say why.
        over = local_client.get(f"{endpoint}?page_size=1001")
        assert over.status_code == 400
        assert "page_size may be at most 1000" in over.json()["error"]["message"]
        past = local_client.get(f"{endpoint}?page=5&page_size=10")
        assert past.status_code == 400
        assert "past the last page (0)" in past.json()["error"]["message"]
        # max boundary accepted, and page 0 of an empty list is no error
        page = local_client.get(f"{endpoint}?page=0&page_size=1000").json()
        assert page["pagination"]["page"] == 0
        assert page["pagination"]["page_size"] == 1000

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
        data = local_client.get(f"{endpoint}?{filter_param}&page=0&page_size=10").json()
        assert data["results"] == []
        assert data["pagination"]["total"] == 0
        assert data["pagination"]["page_size"] == 10


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

        mixed = local_client.get("/v4/genomes?q=test&search_fields=digest,invalid_field,species")
        assert "Invalid search fields: invalid_field" in mixed.json()["error"]["message"]

        multi = local_client.get("/v4/genomes?q=test&search_fields=invalid1,invalid2")
        assert "Invalid search fields: invalid1, invalid2" in multi.json()["error"]["message"]


class TestLocalRouteOwnership:
    """Each path is defined by one handler on one router.

    See docs/refgenie/design-notes.md ("Route ownership across app modes"): the
    alias and asset-files paths are defined once, on the shared router, which
    both modes include at /v4.
    """

    def test_local_mode_serves_alias_endpoints(self, local_client):
        data = local_client.get("/v4/aliases").json()
        assert "results" in data
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

    from refgenie.server.routers import data_channel

    async def fake_fetch(index_url):
        return NATIVE_INDEX

    # The default config ships a 'registry' channel pointing at the Pages URL;
    # stub the network fetch so tests are hermetic but redirects still use the
    # real configured index_url.
    monkeypatch.setattr(data_channel, "_fetch_index", fake_fetch)

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


REGISTRY_INDEX_URL = "https://refgenie.github.io/refgenie-registry/index.yaml"


@pytest.mark.parametrize(
    "path, expected_tail",
    [
        ("recipes/%C3%A9.yaml", "recipes/%C3%A9.yaml"),
        ("recipes/a%20b.yaml", "recipes/a%20b.yaml"),
    ],
)
def test_redirect_reencodes_unsafe_path_characters(data_channel_client, path, expected_tail):
    """Non-latin-1 and space characters are percent-encoded in Location, not a 500."""
    r = data_channel_client.get(f"/v4/data_channel/registry/{path}", follow_redirects=False)
    assert r.status_code == 302
    assert r.headers["location"] == (
        f"https://refgenie.github.io/refgenie-registry/{expected_tail}"
    )


def test_nested_recipe_path_redirects_to_upstream(data_channel_client):
    r = data_channel_client.get(
        "/v4/data_channel/registry/recipes/bwa_index/recipe.yaml",
        follow_redirects=False,
    )
    assert r.status_code == 302
    # urljoin against the configured index_url replaces the trailing index.yaml.
    assert r.headers["location"] == (
        "https://refgenie.github.io/refgenie-registry/recipes/bwa_index/recipe.yaml"
    )
    assert "__" not in r.headers["location"]


def test_asset_class_path_redirects_to_upstream(data_channel_client):
    r = data_channel_client.get(
        "/v4/data_channel/registry/asset_classes/bwa_index.yaml",
        follow_redirects=False,
    )
    assert r.status_code == 302
    assert r.headers["location"] == (
        "https://refgenie.github.io/refgenie-registry/asset_classes/bwa_index.yaml"
    )


def test_unknown_channel_is_404(data_channel_client):
    assert data_channel_client.get("/v4/data_channel/nope/index.yaml").status_code == 404
    assert data_channel_client.get("/v4/data_channel/nope/recipes/x/recipe.yaml").status_code == 404


# --- resolve_channel_file_url: pure-function coverage of the redirect guard --


@pytest.mark.parametrize(
    "file_path, expected_suffix",
    [
        ("recipes/bwa_index/recipe.yaml", "recipes/bwa_index/recipe.yaml"),
        ("asset_classes/bwa_index.yaml", "asset_classes/bwa_index.yaml"),
        ("index.yaml", "index.yaml"),
    ],
)
def test_resolve_channel_file_url_allows_relative_paths(file_path, expected_suffix):
    target = resolve_channel_file_url(REGISTRY_INDEX_URL, file_path)
    assert target == f"https://refgenie.github.io/refgenie-registry/{expected_suffix}"


@pytest.mark.parametrize(
    "file_path",
    [
        "//evil.example/a",
        "///evil.example/a",
        "https://evil.example/a",
        "http://evil.example",
        "https:evil.example/a",
        "javascript:alert(1)",
        "/abs/path",
        "../x",
        "recipes/../../x",
        "recipes/./x",
        "recipes/%2e%2e/%2e%2e/x",
        "%2E%2E/x",
        "\\evil.example/a",
        "recipes\\..\\x",
        "recipes/x?y=1",
        "recipes/x#frag",
        "recipes/x\r\nSet-Cookie: a=b",
        "",
    ],
)
def test_resolve_channel_file_url_rejects_unsafe_paths(file_path):
    assert resolve_channel_file_url(REGISTRY_INDEX_URL, file_path) is None


def test_resolve_channel_file_url_keeps_port():
    base = "https://up.example:8443/ch/index.yaml"
    target = resolve_channel_file_url(base, "a.yaml")
    assert target == "https://up.example:8443/ch/a.yaml"


# --- route-level: the guard applied through the real redirect endpoint ------


@pytest.mark.parametrize(
    "path",
    [
        "/v4/data_channel/registry/%2F%2Fevil.example/a",
        "/v4/data_channel/registry///evil.example/a",
        "/v4/data_channel/registry/https://evil.example/a",
        "/v4/data_channel/registry/https%3A%2F%2Fevil.example%2Fa",
        "/v4/data_channel/registry/%2E%2E/%2E%2E/x",
        "/data_channel/registry/%2F%2Fevil.example/a",
    ],
)
def test_redirect_route_rejects_unsafe_paths(data_channel_client, path):
    r = data_channel_client.get(path, follow_redirects=False)
    assert r.status_code == 404
    assert "location" not in r.headers


# --- DataChannel index_url validation -------------------------------------


def test_data_channel_rejects_non_http_index_url():
    with pytest.raises(ValidationError):
        DataChannel(name="x", protocol="https", index_url="ftp://h/index.yaml")


def test_data_channel_rejects_index_url_without_host():
    with pytest.raises(ValidationError):
        DataChannel(name="x", protocol="https", index_url="index.yaml")


def test_data_channel_rejects_protocol_mismatch():
    with pytest.raises(ValidationError):
        DataChannel(name="x", protocol="http", index_url="https://h/index.yaml")


class TestEndpointHitCounter:
    """The hit counter sees routes added via include_router, with their path params."""

    @pytest.fixture
    def client(self, tmp_path, fixtures_path):
        spec = {**_fasta_spec(stage_modes=("archive",)), "tarball_bytes": b"tarball"}
        tc, digests = _server_client(tmp_path, fixtures_path, [spec])
        with tc as client:
            yield client, digests["asset"]

    def test_archive_download_is_counted(self, client):
        client, digest = client
        assert client.get(f"/v4/archives/{digest}/download").status_code == 200

        collector = client.app.state.endpoint_hits_collector
        stats = collector.get_stats("download_archive")
        assert list(stats) == [frozenset({("asset_digest", digest)})]
        assert stats[frozenset({("asset_digest", digest)})].total_hits == 1

    def test_unmatched_path_is_not_counted(self, client):
        client, _ = client
        assert client.get("/v4/nonexistent").status_code == 404
        assert client.app.state.endpoint_hits_collector.get_stats("download_archive") == {}


# --- genome sequences and comparison ------------------------------------------


@pytest.mark.component
class TestGenomeSequences:
    """``/v4/genomes/{digest}/sequences`` and ``/compare/{other}`` over a real store."""

    @pytest.fixture
    def world(self, tmp_path, fixtures_path):
        """Two real genomes: rCRSd (one sequence) and rCRSd-extra (two)."""
        rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        rgc.database.init(genome_folder=tmp_path / "genomes")
        register_fasta(rgc, fixtures_path)
        digest, _ = rgc.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa", description="", alias_names=["rCRSd"]
        )
        other, _ = rgc.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd-extra.fa",
            description="",
            alias_names=["rCRSd_extra"],
        )
        return rgc, digest, other

    @pytest.fixture
    def client(self, world):
        with make_server_client(world[0]) as client:
            yield client

    def test_sequences_from_real_store(self, world, client):
        _, digest, _ = world
        response = client.get(f"/v4/genomes/{digest}/sequences")
        assert response.status_code == 200
        data = response.json()
        assert data["digest"] == digest
        assert data["n_sequences"] == len(data["sequences"]) == 1
        for row in data["sequences"]:
            assert row["name"]
            assert row["length"] > 0
            assert len(row["sha512t24u"]) == 32
            assert not row["sha512t24u"].startswith("SQ.")
            assert len(row["md5"]) == 32
            assert row["alphabet"]
        assert data["total_length"] == sum(r["length"] for r in data["sequences"])
        assert {"names", "lengths", "sequences"} <= set(data["level1"])
        assert {"namespace": "refgenie", "alias": "rCRSd"} in data["collection_aliases"]

    def test_sequence_aliases_are_listed(self, world, client):
        rgc, digest, _ = world
        before = client.get(f"/v4/genomes/{digest}/sequences").json()
        bare = before["sequences"][0]["sha512t24u"]
        rgc.store_router.default_store.add_sequence_alias("test_ns", "chrM_alias", bare)
        row = client.get(f"/v4/genomes/{digest}/sequences").json()["sequences"][0]
        assert {"namespace": "test_ns", "alias": "chrM_alias"} in row["aliases"]

    def test_sequences_unknown_genome_is_404(self, client):
        response = client.get(f"/v4/genomes/{fake_digest('nope')}/sequences")
        assert response.status_code == 404

    def test_compare_self(self, world, client):
        _, digest, _ = world
        data = client.get(f"/v4/genomes/{digest}/compare/{digest}").json()
        assert data["digests"] == {"a": digest, "b": digest}
        counts = data["array_elements"]
        assert counts["a_and_b_count"]["names"] == counts["a_count"]["names"]
        assert all(data["array_elements"]["a_and_b_same_order"].values())
        assert data["summary"] == "Identical."
        assert data["name_map"] == {}

    def test_compare_two_genomes(self, world, client):
        _, digest, other = world
        data = client.get(f"/v4/genomes/{digest}/compare/{other}").json()
        assert data["digests"] == {"a": digest, "b": other}
        assert set(data["attributes"]) == {"a_only", "b_only", "a_and_b"}
        counts = data["array_elements"]
        assert counts["a_count"]["sequences"] == 1
        assert counts["b_count"]["sequences"] == 2
        assert counts["a_and_b_count"]["sequences"] == 1
        # One sentence and the renames, so clients need not fetch both collections.
        assert data["summary"].startswith("All 1 sequences of")
        assert isinstance(data["name_map"], dict)

    def test_compare_unknown_is_404(self, world, client):
        _, digest, _ = world
        response = client.get(f"/v4/genomes/{digest}/compare/{fake_digest('nope')}")
        assert response.status_code == 404

    def test_the_client_compare_is_the_servers_compare(self, world, tmp_path):
        """``refgenie compare`` of a remote genome asks the server through the client."""
        server_rg, digest, other = world
        client_rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        client_rg.database.init(genome_folder=tmp_path / "client_genomes")
        url = "http://mock-refgenie-server"
        with serve_refgenie(client_rg, server_rg, url):
            remote = client_rg.servers.client(url).genome.compare(digest, other)
        # The server adds its wording; the GA4GH comparison underneath is the same.
        assert remote.pop("summary") and isinstance(remote.pop("name_map"), dict)
        assert remote == server_rg.genome.compare(digest, other)

    def test_sequences_in_local_mode(self, world):
        """The local demo (``refgenie dash``) reads the same endpoint."""
        rgc, digest, _ = world
        with make_local_client(rgc) as client:
            response = client.get(f"/v4/genomes/{digest}/sequences")
            missing = client.get(f"/v4/genomes/{fake_digest('nope')}/sequences")
        assert response.status_code == 200
        assert response.json()["n_sequences"] == 1
        assert missing.status_code == 404
        assert missing.json()["error"]["code"] == "genome_not_found"
