"""Genome sources (``refgenie.managers.sources.genomes``) and the genome CLI on top of them."""

import logging
from contextlib import contextmanager, nullcontext
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pytest

from refgenie.managers.sources.genomes import (
    SERVICE_INFO_PATHS,
    RefgenieServerSource,
    RefgetStoreSource,
    clear_source_cache,
    make_source,
    store_url_from_service_info,
)
from tests.support.factories import make_server_rgc
from tests.support.remotes import MockRemoteSource, assert_mock_remote_source_conforms
from tests.support.server import make_server_app, needs_server

assert_mock_remote_source_conforms()

STORE_URL = "https://refgenie.s3.us-east-1.amazonaws.com/refget-store/jungle/"
API_URL = "https://api.refgenie.org"


@contextmanager
def served_by(app):
    """Answer every ``refgenie.utils.http`` request with the real ``app``.

    The TestClient is not entered, so the lifespan never loads the store. The
    store source is mocked; yields that mock.
    """
    from fastapi.testclient import TestClient

    with (
        patch("refgenie.utils.http.make_client", return_value=nullcontext(TestClient(app))),
        patch("refgenie.managers.sources.genomes.RefgetStoreSource") as store,
    ):
        yield store


def api_app(tmp_path, *, store_url=STORE_URL, web_dist=None):
    """A real server app whose store registry holds ``store_url`` (if any)."""
    rgc = make_server_rgc(tmp_path)
    if store_url:
        rgc.store.add("jungle", store_url)
    return make_server_app(rgc, web_dist=web_dist)


@pytest.fixture
def server_source(tmp_path):
    """A RefgenieServerSource bootstrapped from the real app, its store mocked."""
    with served_by(api_app(tmp_path)) as store:
        source = RefgenieServerSource(API_URL)
    source.store = store.return_value
    return source


@contextmanager
def patch_get(mock=None, **kwargs):
    """Fake the GET behind ``refgenie.utils.http.make_client``, for answers the
    real app never gives (network down, bad JSON, a static store host)."""
    get = mock if mock is not None else MagicMock(**kwargs)
    client = MagicMock()
    client.get = get
    client.__enter__.return_value = client
    with patch("refgenie.utils.http.make_client", return_value=client):
        yield get


def _response(status_code=200, json_data=None, content_type="application/json", json_error=False):
    response = MagicMock()
    response.status_code = status_code
    response.headers = {"content-type": content_type}
    if json_error:
        response.json.side_effect = ValueError("not JSON")
    else:
        response.json.return_value = json_data
    return response


def _service_info(store_url=STORE_URL, enabled=True):
    store = {"enabled": enabled}
    if enabled and store_url:
        store["url"] = store_url
    return {
        "id": "org.refgenie.api",
        "name": "Refgenie",
        "seqcol": {"url": "/seqcol", "refget_store": store},
    }


_VALID_MANIFEST = {
    "version": 1,
    "seqdata_path_template": "sequences/%s2/%s.seq",
    "collections_path_template": "collections/%s.rgsi",
    "sequence_index": "sequences.rgsi",
    "collection_index": "collections.rgci",
    "mode": "Encoded",
    "created_at": "2024-01-01T00:00:00Z",
}


# ---------------------------------------------------------------------------
# Server bootstrap: a server URL becomes its RefgetStore by one service-info fetch
# ---------------------------------------------------------------------------


@needs_server
class TestServiceInfoBootstrap:
    def test_server_resolves_to_its_store(self, tmp_path):
        with served_by(api_app(tmp_path)) as store:
            source = RefgenieServerSource(API_URL)
        assert source.store_url == STORE_URL
        assert source.url == API_URL
        assert store.call_args[0][0] == STORE_URL

    def test_server_without_a_store_is_refused(self, tmp_path):
        with served_by(api_app(tmp_path, store_url=None)):
            with pytest.raises(ConnectionError) as exc:
                RefgenieServerSource(API_URL)
        assert API_URL in str(exc.value)
        assert "RefgetStore" in str(exc.value)

    def test_spa_html_for_rgstore_json_falls_through_to_the_server(self, tmp_path):
        """A SPA answers 200 HTML for every path, rgstore.json included."""
        (tmp_path / "web").mkdir()
        (tmp_path / "web" / "index.html").write_text("<!doctype html><div id='root'></div>")
        with served_by(api_app(tmp_path, web_dist=tmp_path / "web")):
            assert isinstance(make_source(API_URL), RefgenieServerSource)

    def test_metadata_calls_go_to_the_store(self, server_source):
        assert server_source.verify_collection("ABC") is server_source.store.verify_collection(
            "ABC"
        )
        server_source.list_collections(page=2, page_size=7)
        server_source.store.list_collections.assert_called_once_with(page=2, page_size=7)


class TestServiceInfoFailures:
    def test_falls_back_to_the_seqcol_service_info(self):
        def get(url):
            if url.endswith("/seqcol/service-info"):
                return _response(json_data=_service_info())
            return _response(status_code=404)

        with (
            patch_get(side_effect=get),
            patch("refgenie.managers.sources.genomes.RefgetStoreSource"),
        ):
            assert RefgenieServerSource(API_URL).store_url == STORE_URL

    def test_no_service_info_anywhere_raises_listing_both_paths(self):
        with patch_get(return_value=_response(status_code=404)):
            with pytest.raises(ConnectionError) as exc:
                RefgenieServerSource("https://nothing.example.com")
        message = str(exc.value)
        for path in SERVICE_INFO_PATHS:
            assert "https://nothing.example.com" + path in message

    def test_unparseable_service_info_logs_a_warning(self, caplog):
        with caplog.at_level(logging.WARNING):
            with patch_get(return_value=_response(json_error=True)):
                with pytest.raises(ConnectionError):
                    RefgenieServerSource("https://broken.example.com")
        assert any("unparseable JSON" in r.message for r in caplog.records)

    def test_network_failure_logs_a_warning(self, caplog):
        with caplog.at_level(logging.WARNING):
            with patch_get(side_effect=OSError("connection refused")):
                with pytest.raises(ConnectionError):
                    RefgenieServerSource("https://down.example.com")
        assert any("connection refused" in r.message for r in caplog.records)

    def test_make_source_error_names_the_store_probe_and_both_service_infos(self):
        with patch_get(return_value=_response(status_code=404)):
            with pytest.raises(ConnectionError) as exc:
                make_source("https://nothing.example.com")
        message = str(exc.value)
        assert "https://nothing.example.com/rgstore.json" in message
        for path in SERVICE_INFO_PATHS:
            assert "https://nothing.example.com" + path in message


class TestStoreUrlParsing:
    def test_reads_the_nested_shape(self):
        assert store_url_from_service_info(_service_info()) == STORE_URL

    @pytest.mark.parametrize(
        "info",
        [
            None,
            "not a dict",
            {},
            {"seqcol": "not a dict"},
            {"seqcol": {}},
            {"seqcol": {"refget_store": None}},
            {"seqcol": {"refget_store": {"enabled": False, "url": STORE_URL}}},
        ],
    )
    def test_returns_none_for_anything_that_does_not_advertise_a_store(self, info):
        assert store_url_from_service_info(info) is None


# ---------------------------------------------------------------------------
# make_source: its cache and RefgetStore detection
# ---------------------------------------------------------------------------


@needs_server
def test_make_source_caches_by_normalized_url(tmp_path):
    with served_by(api_app(tmp_path)):
        first = make_source(API_URL)
        assert make_source(API_URL + "/") is first
        assert make_source("https://other.example.com") is not first
        clear_source_cache()
        assert make_source(API_URL) is not first


class TestRgstoreDetection:
    """A store is detected by its rgstore.json body, whatever the content type."""

    @pytest.mark.parametrize("content_type", ["application/json", "application/octet-stream"])
    def test_valid_manifest_is_a_refget_store(self, tmp_path, content_type):
        response = _response(json_data=_VALID_MANIFEST, content_type=content_type)
        with (
            patch_get(return_value=response),
            patch("refgenie.managers.sources.genomes.RefgetStoreSource") as MockStore,
        ):
            result = make_source("https://store.example.com", cache_dir=tmp_path)
        MockStore.assert_called_once_with("https://store.example.com", tmp_path)
        assert result is MockStore.return_value

    def test_non_manifest_json_falls_through_to_the_server_bootstrap(self):
        with patch_get(return_value=_response(json_data={"version": 1, "hello": "world"})):
            with pytest.raises(ConnectionError, match="does not advertise a RefgetStore"):
                make_source("https://example.com")


class TestRefgetStoreSourceAliases:
    """RefgetStoreSource pulls the collection alias tables on init, since
    ``open_remote`` alone leaves namespace lookups empty."""

    def test_init_pulls_aliases(self):
        mock_store = MagicMock()
        with patch("gtars.refget.RefgetStore") as MockStore:
            MockStore.open_remote.return_value = mock_store
            RefgetStoreSource("https://example.com/store", Path("/tmp/cache"))
        MockStore.open_remote.assert_called_once()
        mock_store.pull_aliases.assert_called_once()

    def test_resolve_alias_in_namespace_after_pull(self):
        mock_store = MagicMock()
        meta = MagicMock()
        meta.digest = "DIGEST123"
        mock_store.get_collection_metadata_by_alias.return_value = meta
        with patch("gtars.refget.RefgetStore") as MockStore:
            MockStore.open_remote.return_value = mock_store
            src = RefgetStoreSource("https://example.com/store", Path("/tmp/cache"))
        assert src.resolve_alias_in_namespace("base", "name") == "DIGEST123"
        mock_store.get_collection_metadata_by_alias.assert_called_once_with("name", "base")

    def test_init_tolerates_pull_aliases_failure(self):
        """A store with no alias tables must still open (pull failure is non-fatal)."""
        mock_store = MagicMock()
        mock_store.pull_aliases.side_effect = OSError("no alias tables")
        with patch("gtars.refget.RefgetStore") as MockStore:
            MockStore.open_remote.return_value = mock_store
            src = RefgetStoreSource("https://example.com/store", Path("/tmp/cache"))
        assert src.url == "https://example.com/store"


# ---------------------------------------------------------------------------
# Remote genomes: fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def remote_digest():
    return "abcdefghijklmnopqrstuvwxyz123456"


@pytest.fixture
def mock_collection_data():
    return {
        "names": ["chr1", "chr2", "chrM"],
        "lengths": [248956422, 242193529, 16569],
        "sequences": ["seq_digest_1", "seq_digest_2", "seq_digest_3"],
    }


@pytest.fixture
def mock_source(remote_digest, mock_collection_data):
    return MockRemoteSource(
        collections={remote_digest: mock_collection_data},
        store_url="https://store.example.com",
    )


@pytest.fixture
def make_remote_source(remote_digest, mock_collection_data):
    """Build a MockRemoteSource; defaults to one collection at `remote_digest`."""

    def _make(*, digests=None, store_url="https://store.example.com"):
        collections = {d: mock_collection_data for d in (digests or [remote_digest])}
        return MockRemoteSource(collections=collections, store_url=store_url)

    return _make


@pytest.fixture
def init_remote(refgenie_minimal):
    """Register a remote genome; returns (digest, created)."""

    def _init(source, alias, digest=None, description=None):
        return refgenie_minimal.genome.initialize_genome(
            source=source,
            digest=digest,
            description=description if description is not None else alias,
            alias_names=[alias],
        )

    return _init


# ---------------------------------------------------------------------------
# Remote genome initialization and operations
# ---------------------------------------------------------------------------


class TestRemoteGenome:
    """Remote genome initialization and operations via RemoteGenomeSource."""

    def test_init_remote_success(self, refgenie_minimal, remote_digest, mock_source, init_remote):
        """Remote init creates a Genome record with the store URL and alias."""
        digest, created = init_remote(
            mock_source, "test_remote", digest=remote_digest, description="Remote test genome"
        )
        assert created is True
        assert digest == remote_digest

        genome = refgenie_minimal.genome.get(digest)
        assert genome.remote_url == "https://store.example.com"
        assert genome.description == "Remote test genome"
        assert refgenie_minimal.alias.resolve("test_remote") == remote_digest

    def test_init_collection_not_found(self, refgenie_minimal, remote_digest):
        """Remote init raises when the collection is absent on the source."""
        source = MockRemoteSource(collections={})
        with pytest.raises(ValueError, match="not found on source"):
            refgenie_minimal.genome.initialize_genome(
                source=source,
                digest=remote_digest,
                description="Missing genome",
                alias_names=["missing"],
            )

    def test_init_no_store_url(
        self, refgenie_minimal, remote_digest, make_remote_source, init_remote
    ):
        """Remote init without a store_url still succeeds, leaving remote_url unset."""
        source = make_remote_source(store_url=None)
        digest, created = init_remote(
            source, "no_store", digest=remote_digest, description="No store URL genome"
        )
        assert created is True
        assert refgenie_minimal.genome.get(digest).remote_url is None

    def test_init_needs_a_digest(self, refgenie_minimal, mock_source):
        """A source without a digest raises; ``ServerManager.locate`` finds digests."""
        with pytest.raises(ValueError, match="needs a source and a digest"):
            refgenie_minimal.genome.initialize_genome(
                source=mock_source,
                description="No digest",
                alias_names=["hg38"],
            )

    def test_init_missing_source(self, refgenie_minimal, remote_digest):
        """A digest without a source raises."""
        with pytest.raises(ValueError, match="needs a source and a digest"):
            refgenie_minimal.genome.initialize_genome(
                digest=remote_digest,
                description="Missing source",
                alias_names=["no_source"],
            )

    def test_init_mutual_exclusion(
        self, refgenie_minimal, remote_digest, rCRSd_fasta_file_path, mock_source
    ):
        """Cannot specify both fasta_file_path and source/digest."""
        with pytest.raises(ValueError, match="Cannot specify both"):
            refgenie_minimal.genome.initialize_genome(
                fasta_file_path=rCRSd_fasta_file_path,
                source=mock_source,
                digest=remote_digest,
                description="Both sources",
                alias_names=["both"],
            )

    def test_init_no_source_at_all(self, refgenie_minimal):
        """Must specify either fasta or source parameters."""
        with pytest.raises(ValueError, match="Must specify either"):
            refgenie_minimal.genome.initialize_genome(
                description="No source", alias_names=["no_source"]
            )

    def test_init_already_exists(self, refgenie_minimal, remote_digest, mock_source, init_remote):
        """Re-initializing an existing genome raises."""
        init_remote(mock_source, "first_alias", digest=remote_digest, description="First")
        with pytest.raises(ValueError, match="already exists"):
            init_remote(mock_source, "second_alias", digest=remote_digest, description="Second")

    def test_compare_both_remote_same_server(
        self, refgenie_minimal, make_remote_source, init_remote
    ):
        """Two remote genomes on the same store compare their level 2 from that store."""
        digest_a = "aaaaaaaaaaaabbbbbbbbbbbbcccccccc"
        digest_b = "ddddddddddddeeeeeeeeeeeeffffffff"
        source = make_remote_source(digests=[digest_a, digest_b])
        for d, name in ((digest_a, "genome_a"), (digest_b, "genome_b")):
            init_remote(source, name, digest=d)

        level2 = {"names": ["chr1"], "lengths": [8], "sequences": ["SQ.x"]}
        with patch.object(
            refgenie_minimal.genome.__class__, "_get_remote_level2", return_value=level2
        ) as remote:
            result = refgenie_minimal.genome.compare(digest_a, digest_b)
        assert {c.args for c in remote.call_args_list} == {
            ("https://store.example.com", digest_a),
            ("https://store.example.com", digest_b),
        }
        assert result["digests"] == {"a": digest_a, "b": digest_b}
        assert result["attributes"]["a_and_b"] == ["lengths", "names", "sequences"]

    def test_compare_one_local_one_remote(
        self,
        refgenie_minimal,
        rCRSd_fasta_file_path,
        mock_collection_data,
        make_remote_source,
        init_remote,
    ):
        """Comparing a local and a remote genome uses local comparison."""
        local_digest, _ = refgenie_minimal.genome.initialize_genome(
            fasta_file_path=rCRSd_fasta_file_path,
            description="Local genome",
            alias_names=["local_genome"],
        )
        remote_digest = "xxxxxxxxyyyyyyyyzzzzzzzzaabbccdd"
        source = make_remote_source(digests=[remote_digest])
        init_remote(source, "remote_comp", digest=remote_digest, description="Remote genome")
        with patch.object(
            refgenie_minimal.genome.__class__,
            "_get_remote_level2",
            return_value=mock_collection_data,
        ):
            result = refgenie_minimal.genome.compare(local_digest, remote_digest)
        assert set(result) >= {"digests", "attributes", "array_elements"}

    def test_list_shows_remote_status(
        self, refgenie_minimal, rCRSd_fasta_file_path, make_remote_source, init_remote
    ):
        """Listing reports remote_url for remote genomes and None for local ones."""
        refgenie_minimal.genome.initialize_genome(
            fasta_file_path=rCRSd_fasta_file_path,
            description="Local",
            alias_names=["local_list"],
        )
        remote_digest = "aabbccdd11223344aabbccdd11223344"
        source = make_remote_source(digests=[remote_digest])
        init_remote(source, "remote_list", digest=remote_digest, description="Remote")
        genome_map = {g.digest: g for g in refgenie_minimal.genome.list_all()}
        local_digest = refgenie_minimal.alias.resolve("local_list")
        assert genome_map[local_digest].remote_url is None
        assert genome_map[remote_digest].remote_url == "https://store.example.com"

    @pytest.mark.parametrize(
        "store_url, expected_calls",
        [
            ("https://store.example.com", "one"),
            (None, "none"),
        ],
        ids=["with-store-url", "without-store-url"],
    )
    def test_import_remote_collection_gated_on_store_url(
        self, refgenie_minimal, remote_digest, make_remote_source, store_url, expected_calls
    ):
        """Remote init imports the collection only when the source exposes a store URL."""
        source = make_remote_source(store_url=store_url)
        with patch.object(
            refgenie_minimal.genome.__class__, "_import_remote_collection"
        ) as mock_import:
            refgenie_minimal.genome.initialize_genome(
                source=source,
                digest=remote_digest,
                description="Import test",
                alias_names=["import_test"],
            )
        expected = [call(store_url, remote_digest)] if store_url else []
        assert mock_import.call_args_list == expected

    def test_import_failure_does_not_block_init(
        self, refgenie_minimal, remote_digest, mock_source, init_remote
    ):
        """A failing collection import does not prevent genome creation."""
        with patch("gtars.refget.RefgetStore") as MockGtars:
            MockGtars.open_remote.side_effect = Exception("Connection failed")
            digest, created = init_remote(
                mock_source, "fail_import", digest=remote_digest, description="Fail import"
            )
        assert created is True
        assert refgenie_minimal.genome.exists(remote_digest)

    def test_list_shows_aliases_for_remote(
        self, refgenie_minimal, remote_digest, mock_source, init_remote
    ):
        """genome list shows aliases for remote-sourced genomes."""
        with patch.object(refgenie_minimal.genome.__class__, "_import_remote_collection"):
            init_remote(
                mock_source,
                "list_alias_test",
                digest=remote_digest,
                description="List alias test",
            )
        assert "list_alias_test" in refgenie_minimal.alias.get_for_genome(remote_digest)

    def test_compare_level2_helper(self):
        """compare_level2 returns the GA4GH seqcol comparison object."""
        from refgenie.utils.seqcol import compare_level2 as _compare_level2

        level2_a = {"names": ["chr1", "chr2"], "lengths": [100, 200], "sequences": ["s1", "s2"]}
        level2_b = {"names": ["chr1", "chr3"], "lengths": [100, 300], "sequences": ["s1", "s3"]}
        result = _compare_level2(level2_a, level2_b, "digest_a", "digest_b")

        assert result["digests"] == {"a": "digest_a", "b": "digest_b"}
        assert sorted(result["attributes"]["a_and_b"]) == ["lengths", "names", "sequences"]
        elements = result["array_elements"]
        assert set(elements) == {"a_count", "b_count", "a_and_b_count", "a_and_b_same_order"}
        assert elements["a_count"]["names"] == 2
        assert elements["b_count"]["names"] == 2
        assert elements["a_and_b_count"]["names"] == 1  # chr1


# ---------------------------------------------------------------------------
# Genome CLI models and browse/sync handlers
# ---------------------------------------------------------------------------


def _browse_cmd(server_url=None, page=0, page_size=20):
    from types import SimpleNamespace

    return SimpleNamespace(server_url=server_url, page=page, page_size=page_size)


def _sync_cmd(server_url=None, page_size=1000):
    from types import SimpleNamespace

    return SimpleNamespace(server_url=server_url, page_size=page_size)


class _Meta:
    """Stand-in for gtars' SequenceCollectionMetadata."""

    def __init__(self, digest, n_sequences, names_digest):
        self.digest = digest
        self.n_sequences = n_sequences
        self.names_digest = names_digest


class TestGenomeCliModels:
    """Validation of the genome CLI pydantic models."""

    def test_sync_model(self):
        """GenomeSyncModel exposes server-url/page-size with correct defaults."""
        from refgenie.cli.commands.genome import GenomeSyncModel

        m = GenomeSyncModel()
        assert m.server_url is None
        assert m.page_size == 1000
        assert not hasattr(m, "prefix")

        m2 = GenomeSyncModel.model_validate({"server-url": "https://x", "page-size": 5})
        assert m2.server_url == "https://x"
        assert m2.page_size == 5


class TestGenomeBrowse:
    """genome browse output rendering."""

    def test_browse_prints_real_counts_not_question_marks(self, refgenie_minimal, capsys):
        """browse must print real digest / n_sequences / names_digest.

        Regression: refget's REST list endpoint flattens results to bare
        digest strings, so ``getattr(str, "n_sequences", "?")`` printed
        "?  ? sequences" for every row. Store-backed sources return objects
        carrying all three fields.
        """
        from refgenie.cli.commands import genome as handlers

        source = MagicMock()
        source.list_collections.return_value = {
            "results": [_Meta("DIGEST_A", 25, "NAMES_A"), _Meta("DIGEST_B", 3, "NAMES_B")]
        }
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_browse(_browse_cmd(server_url="https://s"), refgenie_minimal)

        out = capsys.readouterr().out
        assert "?" not in out
        assert "DIGEST_A\t25 sequences\tNAMES_A" in out
        assert "DIGEST_B\t3 sequences\tNAMES_B" in out

    def test_store_list_collections_supplies_browse_fields(
        self, tmp_path, refgenie_minimal, capsys
    ):
        """The real gtars store returns objects with the three fields browse reads."""
        from gtars.refget import RefgetStore

        from refgenie.cli.commands import genome as handlers

        fasta = tmp_path / "tiny.fa"
        fasta.write_text(">chr1\nACGTACGTAC\n>chr2\nTTTTGGGGCC\n")

        store = RefgetStore.in_memory()
        store.add_sequence_collection_from_fasta(str(fasta))

        source = MagicMock()
        source.list_collections.return_value = store.list_collections()
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_browse(_browse_cmd(server_url="https://s"), refgenie_minimal)

        out = capsys.readouterr().out
        assert "?" not in out
        assert "2 sequences" in out
        digest, n_seqs, names_digest = out.strip().split("\t")
        assert len(digest) == 32
        assert n_seqs == "2 sequences"
        assert len(names_digest) == 32


class TestGenomeSync:
    """genome sync bulk-registration behavior."""

    def test_registers_all_collections(self, refgenie_minimal, mock_collection_data):
        """sync bulk-registers every collection advertised by the source."""
        from refgenie.cli.commands import genome as handlers

        digests = [f"sync{i:028d}" for i in range(3)]
        source = MockRemoteSource(
            collections={d: mock_collection_data for d in digests},
            store_url="https://store.example.com",
        )
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(
                _sync_cmd(server_url="https://server.example.com"), refgenie_minimal
            )
        for d in digests:
            assert refgenie_minimal.genome.get(d).digest == d

    def test_is_idempotent(self, refgenie_minimal, mock_collection_data):
        """Running sync twice does not error and keeps every genome."""
        from refgenie.cli.commands import genome as handlers

        digests = [f"sync{i:028d}" for i in range(2)]
        source = MockRemoteSource(
            collections={d: mock_collection_data for d in digests},
            store_url="https://store.example.com",
        )
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)
        for d in digests:
            assert refgenie_minimal.genome.get(d).digest == d

    def test_falls_back_to_subscriptions(self, refgenie_minimal, mock_collection_data):
        """With no explicit URL, sync uses the configured server subscriptions."""
        from refgenie.cli.commands import genome as handlers

        digests = [f"subs{i:028d}" for i in range(2)]
        source = MockRemoteSource(
            collections={d: mock_collection_data for d in digests},
            store_url="https://store.example.com",
        )
        with (
            patch.object(refgenie_minimal.servers, "subscriptions", return_value=["https://sub"]),
            patch.object(handlers, "make_source", return_value=source),
        ):
            handlers.handle_genome_sync(_sync_cmd(server_url=None), refgenie_minimal)
        for d in digests:
            assert refgenie_minimal.genome.get(d).digest == d

    def test_no_url_no_subscriptions_exits(self, refgenie_minimal):
        """With no URL and no subscriptions, sync exits with an error."""
        from refgenie.cli.commands import genome as handlers

        with patch.object(refgenie_minimal.servers, "subscriptions", return_value=[]):
            with pytest.raises(SystemExit):
                handlers.handle_genome_sync(_sync_cmd(server_url=None), refgenie_minimal)

    def test_carries_name_namespace_aliases_only(self, refgenie_minimal, mock_collection_data):
        """Synced genomes carry the store's 'name' aliases; other namespaces stay out."""
        from refgenie.cli.commands import genome as handlers

        digest = "aliascarry0000000000000000000000"
        source = MockRemoteSource(
            collections={digest: mock_collection_data},
            collection_aliases={
                digest: [
                    ("name", "hg38-broad"),
                    ("genome_assembly", "hg38"),
                    ("accession", "GCA_000001405.15"),
                ]
            },
        )
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)

        assert refgenie_minimal.alias.get_for_genome(digest) == ["hg38-broad"]
        assert not refgenie_minimal.alias.exists("hg38")
        assert not refgenie_minimal.alias.exists("GCA_000001405.15")

    def test_skips_names_that_are_not_valid_aliases(self, refgenie_minimal, mock_collection_data):
        """An iGenomes path such as 'Homo_sapiens/NCBI/GRCh38' is not an alias here."""
        from refgenie.cli.commands import genome as handlers

        digest = "aliasslash000000000000000000000x"
        source = MockRemoteSource(
            collections={digest: mock_collection_data},
            collection_aliases={digest: [("name", "Homo_sapiens/NCBI/GRCh38"), ("name", "hg38")]},
        )
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)

        assert refgenie_minimal.alias.get_for_genome(digest) == ["hg38"]

    def test_never_repoints_existing_alias(
        self, refgenie_minimal, rCRSd_fasta_file_path, mock_collection_data
    ):
        """An alias owned by another genome is skipped, never stolen."""
        from refgenie.cli.commands import genome as handlers

        built_digest, _ = refgenie_minimal.genome.initialize_genome(
            fasta_file_path=rCRSd_fasta_file_path,
            description="Built",
            alias_names=["rCRSd"],
        )
        synced_digest = "aliassteal0000000000000000000000"
        source = MockRemoteSource(
            collections={synced_digest: mock_collection_data},
            collection_aliases={synced_digest: [("name", "rCRSd"), ("name", "other-name")]},
        )
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)

        assert refgenie_minimal.alias.resolve("rCRSd") == built_digest
        assert refgenie_minimal.alias.resolve("other-name") == synced_digest

    def test_applies_fhr_metadata(
        self, refgenie_minimal, rCRSd_fasta_file_path, mock_collection_data
    ):
        """A collection with an FHR sidecar gets the faceted metadata columns.

        The collection is ingested into the local store first (as the real
        sync's collection import does) because apply_fhr also writes the store
        FHR sidecar, which requires the collection to exist there.
        """
        from refgenie.cli.commands import genome as handlers

        metadata, _ = refgenie_minimal.refget_store.add_sequence_collection_from_fasta(
            rCRSd_fasta_file_path
        )
        digest = metadata.digest
        source = MockRemoteSource(
            collections={digest: mock_collection_data},
            fhr={
                digest: {
                    "genome": "hg38",
                    "taxon": {
                        "name": "Homo sapiens",
                        "uri": "https://identifiers.org/taxonomy:9606",
                    },
                    "documentation": "GRCh38 from Broad",
                    "assemblySource": "broad",
                    "accessionID": {"name": "GCA_000001405.15"},
                }
            },
        )
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)

        genome = refgenie_minimal.genome.get(digest)
        assert genome.species_name == "Homo sapiens"
        assert genome.taxon_id == 9606
        assert genome.assembly_source == "broad"
        assert genome.assembly_accession == "GCA_000001405.15"
        assert genome.description == "GRCh38 from Broad"

    def test_no_fhr_leaves_columns_null(self, refgenie_minimal, mock_collection_data):
        """A collection with no FHR still registers, with metadata columns null."""
        from refgenie.cli.commands import genome as handlers

        digest = "fhrnull0000000000000000000000000"
        source = MockRemoteSource(collections={digest: mock_collection_data})
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)

        genome = refgenie_minimal.genome.get(digest)
        assert genome.species_name is None
        assert genome.taxon_id is None

    def test_resync_attaches_new_aliases(self, refgenie_minimal, mock_collection_data):
        """A second sync attaches aliases that appeared at the source since the first."""
        from refgenie.cli.commands import genome as handlers

        digest = "aliaslater0000000000000000000000"
        source = MockRemoteSource(collections={digest: mock_collection_data})
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)
        assert refgenie_minimal.alias.get_for_genome(digest) == []

        source._collection_aliases = {digest: [("name", "late-alias")]}
        with patch.object(handlers, "make_source", return_value=source):
            handlers.handle_genome_sync(_sync_cmd(server_url="https://s"), refgenie_minimal)
        assert refgenie_minimal.alias.get_for_genome(digest) == ["late-alias"]
