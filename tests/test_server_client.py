"""Unit tests for the client half of refgenie's remote features.

Covers LocalMode/ServerMode store selection and RefgenieserverClient
pagination/construction, plus the download progress bar.

NOTE: ``TestDownloadWithProgress`` below is the ONE sanctioned hand-rolled
FastAPI app in the suite (see CLAUDE.md, "Server Endpoint Testing").
Everywhere else, tests must wrap a real Refgenie instance in ``create_app``
(helpers `make_server_app(rg)` / `serve_refgenie(...)` in tests/helpers.py).
The exemption exists because the real ``FileResponse`` archive route cannot
produce a no-``Content-Length`` response, which that client progress-bar
regression test requires.

The push / catalog-transfer (component) half lives in tests/test_remote_push.py.
"""

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import httpx
import pytest

from refgenie import Refgenie
from refgenie.managers.sources.client import RefgenieserverClient
from tests.helpers import real_client_no_init


# ===========================================================================
# LocalMode vs ServerMode store selection
# ===========================================================================


class TestRemoteStore:
    """Federated store router wiring (local on-disk vs server remote)."""

    def test_router_opens_remote_store(self, engine):
        """An enabled remote Store row opens via RefgetStore.open_remote."""
        from refgenie.db.tables import StoreType

        r = Refgenie(database_engine=engine, suppress_migrations=True, server_mode=True)
        r._create_db_and_tables()
        r.store.add("s1", "https://store.example.com", StoreType.remote, priority=10)

        mock_store = MagicMock()
        with patch("refgenie.core.store_router.RefgetStore") as MockRefgetStore:
            MockRefgetStore.open_remote.return_value = mock_store
            store = r.refget_store

        MockRefgetStore.open_remote.assert_called_once()
        assert str(MockRefgetStore.open_remote.call_args[0][1]) == "https://store.example.com"
        assert store is mock_store

    def test_router_opens_on_disk_store_local_mode(self, engine, tmp_path):
        """Local mode's single store opens via RefgetStore.on_disk."""
        r = Refgenie(database_engine=engine, suppress_migrations=True)
        r.init(genome_folder=tmp_path / "genomes")
        mock_store = MagicMock()
        with patch("refgenie.core.store_router.RefgetStore") as MockRefgetStore:
            MockRefgetStore.on_disk.return_value = mock_store
            store = r.refget_store

        MockRefgetStore.on_disk.assert_called_once()
        assert store is mock_store

    def test_getseq_server_mode_raises_on_missing_sequence(self, engine, tmp_path):
        """In server mode, getseq raises RuntimeError when the sequence is missing."""
        from refgenie.db.tables import Alias, Genome

        r = Refgenie(
            database_engine=engine, suppress_migrations=True, server_mode=True
        )
        r.init(genome_folder=tmp_path / "genomes")
        with r._database_session as session:
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

            def get_store(self, name):
                return self.default_store

        r._store_router = _FakeRouter(mock_store)

        with pytest.raises(RuntimeError, match="not available in server mode"):
            r.getseq("test_genome", "chr1")


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
                {"items": [{"id": 1}, {"id": 2}]},
                httpx.ReadTimeout("page 2 timed out"),
            ]
        )
        with pytest.raises(httpx.ReadTimeout):
            client.get_paginated("some_operation", page_size=2)

    def test_partial_pagination_non_strict_returns_partial(self):
        """strict=False keeps the old lenient behaviour, but logs a warning."""
        client = real_client_no_init()
        client.get = MagicMock(
            side_effect=[
                {"items": [{"id": 1}, {"id": 2}]},
                httpx.ReadTimeout("page 2 timed out"),
            ]
        )
        results = client.get_paginated("some_operation", page_size=2, strict=False)
        assert len(results) == 2

    def test_complete_pagination_returns_all_pages(self):
        """A clean run still walks every page."""
        client = real_client_no_init()
        client.get = MagicMock(
            side_effect=[{"items": [{"id": 1}, {"id": 2}]}, {"items": [{"id": 3}]}]
        )
        results = client.get_paginated("some_operation", page_size=2)
        assert [r["id"] for r in results] == [1, 2, 3]

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


# ===========================================================================
# Download progress bar handling (client behavior)
# ===========================================================================


class TestDownloadWithProgress:
    """download_with_progress with and without a Content-Length header.

    NOTE: this is the ONE sanctioned hand-rolled FastAPI app in the suite.
    Everywhere else, tests must wrap a real Refgenie instance in
    ``create_app`` (see tests/helpers.py::serve_refgenie). It survives here
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
