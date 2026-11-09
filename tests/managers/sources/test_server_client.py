"""The client side of remote features: store-mode selection, RefgenieserverClient
pagination and construction, and the download progress bar.

``TestDownloadWithProgress`` is the one hand-rolled FastAPI app in the suite:
the real archive route always sends ``Content-Length``, and that test needs a
response without one.
"""

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import httpx
import pytest

from refgenie import Refgenie
from refgenie.managers.sources.api_ids import (
    API_ID_ALIAS_DIGEST,
    API_ID_ASSET_FILE_DOWNLOAD,
    API_ID_GENOME_ATTRS,
)
from refgenie.managers.sources.client import AliasClient, GenomeClient, RefgenieserverClient
from tests.support.remotes import real_client_no_init
from tests.support.server import make_server_app, requires_server, server_with_asset


# ===========================================================================
# LocalMode vs ServerMode store selection
# ===========================================================================


class TestRemoteStore:
    """Federated store router wiring (local on-disk vs server remote)."""

    def test_router_opens_remote_store(self, engine):
        """Server mode opens a remote Store row, loads it fully, and freezes it."""
        from refgenie.db.tables import StoreType

        r = Refgenie(database_engine=engine, suppress_migrations=True, server_mode=True)
        r.database.create_tables()
        r.store.add("s1", "https://store.example.com", StoreType.remote, priority=10)

        mock_store = MagicMock()
        with patch("refgenie.core.store_router.RefgetStore") as MockRefgetStore:
            MockRefgetStore.open_remote.return_value = mock_store
            store = r.refget_store

        MockRefgetStore.open_remote.assert_called_once()
        assert str(MockRefgetStore.open_remote.call_args[0][1]) == "https://store.example.com"
        # Server mode loads everything once and serves the frozen copy.
        mock_store.load_all_collections.assert_called_once_with()
        mock_store.pull_fhr.assert_called_once_with()
        assert store is mock_store.into_readonly.return_value

    def test_router_opens_on_disk_store_local_mode(self, engine, tmp_path):
        """Local mode's single store opens via RefgetStore.on_disk."""
        r = Refgenie(database_engine=engine, suppress_migrations=True)
        r.database.init(genome_folder=tmp_path / "genomes")
        mock_store = MagicMock()
        with patch("refgenie.core.store_router.RefgetStore") as MockRefgetStore:
            MockRefgetStore.on_disk.return_value = mock_store
            store = r.refget_store

        MockRefgetStore.on_disk.assert_called_once()
        assert store is mock_store

    def test_getseq_server_mode_raises_on_missing_sequence(self, engine, tmp_path):
        """In server mode, getseq raises RuntimeError when the sequence is missing."""
        from refgenie.db.tables import Alias, Genome

        r = Refgenie(database_engine=engine, suppress_migrations=True, server_mode=True)
        r.database.init(genome_folder=tmp_path / "genomes")
        with r.database._database_session as session:
            genome = Genome(digest="testdigest1234567890abcdef123456", description="test")
            session.add(genome)
            session.add(Alias(name="test_genome", genome_digest=genome.digest))
            session.commit()

        mock_store = MagicMock()
        mock_store.get_sequence_by_name.side_effect = KeyError("not found")

        class _FakeRouter:
            names: list = []

            def __init__(self, store):
                self.default_store = store

            def store_for_genome(self, genome_digest, store_name):
                return self.default_store

        r._store_router = _FakeRouter(mock_store)

        with pytest.raises(RuntimeError, match="not available in server mode"):
            r.sequence.get("test_genome", "chr1")


# ===========================================================================
# RefgenieserverClient pagination and construction
# ===========================================================================


class TestServerClientPagination:
    """RefgenieserverClient.get_paginated and construction error handling."""

    def test_partial_pagination_raises(self):
        """A mid-pagination failure must surface, not silently truncate.

        Regression: get_paginated swallowed the exception and returned the pages
        it had, so callers could not tell 2 items from 2000.
        """
        client = real_client_no_init()
        client.get = MagicMock(
            side_effect=[
                {"results": [{"id": 1}, {"id": 2}], "pagination": {"total": 3}},
                httpx.ReadTimeout("page 2 timed out"),
            ]
        )
        with pytest.raises(httpx.ReadTimeout):
            client.get_paginated("some_operation", page_size=2)

    def test_complete_pagination_returns_all_pages(self):
        """A clean run still walks every page."""
        client = real_client_no_init()
        client.get = MagicMock(
            side_effect=[
                {"results": [{"id": 1}, {"id": 2}], "pagination": {"total": 3}},
                {"results": [{"id": 3}], "pagination": {"total": 3}},
            ]
        )
        results = client.get_paginated("some_operation", page_size=2)
        assert [r["id"] for r in results] == [1, 2, 3]
        pages = [call.kwargs["params"] for call in client.get.call_args_list]
        assert pages == [{"page": 0, "page_size": 2}, {"page": 1, "page_size": 2}]

    def test_pagination_stops_at_the_total(self):
        """A full last page is the end: no request for a page past it, which is a 400."""
        client = real_client_no_init()
        client.get = MagicMock(
            return_value={"results": [{"id": 1}, {"id": 2}], "pagination": {"total": 2}}
        )
        assert len(client.get_paginated("some_operation", page_size=2)) == 2
        assert client.get.call_count == 1

    def test_unreachable_server_raises_at_construction(self):
        """An unreachable server must not yield a silently broken client.

        Regression: __init__ caught httpx.ConnectError, logged, and handed back
        a fully-constructed object that failed later with an unrelated message.
        """
        transport = httpx.MockTransport(
            lambda request: (_ for _ in ()).throw(httpx.ConnectError("no route"))
        )
        http_client = httpx.Client(transport=transport)
        with pytest.raises(httpx.ConnectError):
            RefgenieserverClient("http://unreachable.example.com", http_client=http_client)


class TestResourceClientLookups:
    """``alias.resolve`` and ``genome.get`` read one record each."""

    def test_alias_resolve_returns_the_digest(self):
        client = real_client_no_init()
        client.get = MagicMock(return_value={"digest": "abc"})
        assert AliasClient(client).resolve("hg38") == "abc"
        client.get.assert_called_once_with(
            operation_id=API_ID_ALIAS_DIGEST, url_format_params={"name": "hg38"}
        )

    def test_alias_resolve_without_a_digest_is_none(self):
        client = real_client_no_init()
        client.get = MagicMock(return_value="not json")
        assert AliasClient(client).resolve("hg38") is None

    def test_genome_get_returns_the_record(self):
        client = real_client_no_init()
        client.get = MagicMock(return_value={"description": "hg"})
        assert GenomeClient(client).get("abc") == {"description": "hg"}
        client.get.assert_called_once_with(
            operation_id=API_ID_GENOME_ATTRS, url_format_params={"digest": "abc"}
        )

    def test_genome_get_non_dict_is_empty(self):
        client = real_client_no_init()
        client.get = MagicMock(return_value="oops")
        assert GenomeClient(client).get("abc") == {}


@pytest.mark.component
class TestServerClientGet:
    """RefgenieserverClient.get against the real server app."""

    def test_text_response_comes_back_as_text(self, tmp_path, fixtures_path):
        """A text response comes back as text, not JSON-parsed."""
        requires_server()
        from fastapi.testclient import TestClient

        rgc, digests = server_with_asset(
            tmp_path,
            fixtures_path=fixtures_path,
            assets=[
                {
                    "asset_group_name": "fasta",
                    "asset_class": "fasta",
                    "files": {"test.fa": ">chr1\nACGT\n"},
                    "stage_modes": ("file",),
                }
            ],
        )
        with TestClient(make_server_app(rgc)) as tc:
            client = RefgenieserverClient("http://testserver", http_client=tc)
            served = client.get(
                API_ID_ASSET_FILE_DOWNLOAD,
                url_format_params={"asset_digest": digests["fasta"], "file_path": "test.fa"},
            )
        assert served == ">chr1\nACGT\n"


# ===========================================================================
# Download progress bar handling (client behavior)
# ===========================================================================


class TestDownloadWithProgress:
    """download_with_progress with and without a Content-Length header.

    NOTE: this is the ONE sanctioned hand-rolled FastAPI app in the suite.
    Everywhere else, tests must wrap a real Refgenie instance in
    ``create_app`` (see tests/support/server.py::serve_refgenie). It survives here
    because the no-Content-Length case cannot be produced by the real archive
    route (FileResponse always sets Content-Length), and this is a CLIENT
    progress-bar regression test (total=None), not a server-route test.
    """

    @pytest.fixture
    def download_client(self, tmp_path):
        """RefgenieserverClient with injected TestClient for a minimal download app."""
        pytest.importorskip("fastapi", reason="fastapi not installed (dash extras)")
        from fastapi import FastAPI
        from fastapi.responses import Response, StreamingResponse
        from fastapi.testclient import TestClient

        app = FastAPI(openapi_url=None)

        @app.get("/openapi.json")
        def openapi():
            return {
                "openapi": "3.0.0",
                "info": {"title": "Test", "version": "1.0.0", "description": "Test"},
                "tags": [],
                "paths": {
                    "/download/{id}": {"get": {"operationId": "download_op"}},
                    "/download-stream/{id}": {"get": {"operationId": "download_stream_op"}},
                    "/download-fail/{id}": {"get": {"operationId": "download_fail_op"}},
                },
            }

        @app.get("/download/{id}")
        def download(id: str):
            content = b"test content " * 100
            return Response(content=content, headers={"Content-Length": str(len(content))})

        @app.get("/download-stream/{id}")
        def download_stream(id: str):
            return StreamingResponse((b"test content " for _ in range(100)))

        @app.get("/download-fail/{id}")
        def download_fail(id: str):
            """A response that starts, sends some bytes, then blows up mid-stream --
            simulating an interrupted/failed download after headers are sent."""

            def gen():
                yield b"partial content"
                raise RuntimeError("simulated mid-stream failure")

            return StreamingResponse(gen())

        with TestClient(app) as tc:
            yield RefgenieserverClient("http://testserver", http_client=tc)

    def test_download_progress(self, download_client, tmp_path):
        """Download with and without Content-Length works correctly."""
        # -- With Content-Length --
        path1 = tmp_path / "out1.bin"
        download_client.download_with_progress(
            "download_op", path1, url_format_params={"id": "x"}, name="test"
        )
        assert path1.exists() and b"test content" in path1.read_bytes()

        # -- Without Content-Length (bug fix regression: total=None) --
        path2 = tmp_path / "out2.bin"
        download_client.download_with_progress(
            "download_stream_op", path2, url_format_params={"id": "x"}, name="test"
        )
        assert path2.exists() and b"test content" in path2.read_bytes()

    def test_download_renames_from_temp_file_on_success(self, download_client, tmp_path):
        """A successful download lands via an atomic rename, not a direct write.

        Regression: output_path.open("wb") wrote straight to the final
        destination, so an interrupted download could leave a truncated file
        sitting where a complete one belongs.
        """
        path = tmp_path / "out.bin"
        tmp_download_path = path.with_name(path.name + ".part")
        seen_tmp_path_at_rename = {}
        real_replace = os.replace

        def spy_replace(src, dst):
            # At the moment of the rename, the download must have landed at
            # the temp path (not the final path), fully written.
            seen_tmp_path_at_rename["src_exists"] = Path(src).exists()
            seen_tmp_path_at_rename["dst_missing"] = not Path(dst).exists()
            real_replace(src, dst)

        with patch("refgenie.managers.sources.client.os.replace", side_effect=spy_replace):
            download_client.download_with_progress(
                "download_op", path, url_format_params={"id": "x"}, name="test"
            )

        assert seen_tmp_path_at_rename == {"src_exists": True, "dst_missing": True}
        assert path.exists() and b"test content" in path.read_bytes()
        assert not tmp_download_path.exists()

    def test_failed_download_leaves_no_partial_file(self, download_client, tmp_path):
        """A download that fails mid-stream must not leave a truncated file at
        the final path, nor a stray temp file behind.
        """
        path = tmp_path / "out.bin"
        tmp_download_path = path.with_name(path.name + ".part")

        with pytest.raises(Exception):
            download_client.download_with_progress(
                "download_fail_op", path, url_format_params={"id": "x"}, name="test"
            )

        assert not path.exists()
        assert not tmp_download_path.exists()
