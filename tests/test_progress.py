"""The library progress hook, ``refgenie.progress``.

With no sink installed ``emit`` must cost nothing; with one installed the sink
is thread-scoped, and ``download_with_progress`` reports to it instead of
building a rich display.
"""

import threading

import pytest

from tests.helpers import requires_server

requires_server()

from fastapi import FastAPI  # noqa: E402  (must follow the extras guard)
from fastapi.testclient import TestClient  # noqa: E402

from refgenie import progress  # noqa: E402
from refgenie.progress import ProgressAborted  # noqa: E402


class TestProgressHook:
    """With no sink installed ``emit`` must cost nothing (every CLI call); with
    one installed the sink is thread-scoped and always uninstalled on exit,
    because pooled worker threads would otherwise cross-contaminate jobs."""

    def test_emit_without_sink_is_a_noop(self):
        assert progress.active_sink() is None
        assert progress.emit("progress", "x", current=1, total=2) is None

    def test_active_sink_reports_the_installed_sink(self):
        def sink(event):
            pass

        with progress.use_sink(sink):
            assert progress.active_sink() is sink
        assert progress.active_sink() is None

    def test_events_carry_every_field(self):
        events = []
        with progress.use_sink(events.append):
            progress.emit("progress", "downloading", current=5, total=10, unit="bytes", x=1)
        assert len(events) == 1
        event = events[0]
        assert (event.type, event.message, event.current, event.total, event.unit) == (
            "progress",
            "downloading",
            5,
            10,
            "bytes",
        )
        assert event.extra == {"x": 1}

    def test_two_threads_see_only_their_own_sink(self):
        """A ContextVar starts empty in each new thread -- exactly the per-job
        isolation the job manager relies on."""
        results = {"a": [], "b": []}
        started = threading.Barrier(2)

        def worker(name):
            with progress.use_sink(lambda e: results[name].append(e.message)):
                started.wait(timeout=5)
                for i in range(5):
                    progress.emit("progress", f"{name}{i}")

        threads = [threading.Thread(target=worker, args=(n,)) for n in ("a", "b")]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5)
        assert results["a"] == [f"a{i}" for i in range(5)]
        assert results["b"] == [f"b{i}" for i in range(5)]

    def test_token_is_reset_even_on_exception(self):
        with pytest.raises(RuntimeError):
            with progress.use_sink(lambda e: None):
                raise RuntimeError("boom")
        assert progress.active_sink() is None

    def test_sink_exceptions_are_not_swallowed(self):
        """``ProgressAborted`` propagating out of ``emit`` IS the cancellation
        mechanism -- swallowing sink errors would silently disable it."""

        def sink(event):
            raise ProgressAborted("stop")

        with progress.use_sink(sink):
            with pytest.raises(ProgressAborted):
                progress.emit("progress", current=1)


class TestDownloadWithProgressSink:
    """``download_with_progress`` reports to a sink instead of building a rich
    bar. Not constructing the rich display is what makes concurrent pulls
    possible: rich permits one live display per console."""

    @pytest.fixture
    def download_client(self, tmp_path):
        from fastapi.responses import Response

        from refgenie.managers.sources.client import RefgenieserverClient

        app = FastAPI(openapi_url=None)

        @app.get("/openapi.json")
        def openapi():
            return {
                "openapi": "3.0.0",
                "info": {"title": "Test", "version": "1.0.0", "description": "Test"},
                "tags": [],
                "paths": {"/download/{id}": {"get": {"operationId": "download_op"}}},
            }

        @app.get("/download/{id}")
        def download(id: str):
            content = b"x" * 5000
            return Response(content=content, headers={"Content-Length": str(len(content))})

        with TestClient(app) as tc:
            yield RefgenieserverClient("http://testserver", http_client=tc)

    def test_sink_receives_byte_progress_and_no_rich_display(
        self, download_client, tmp_path, monkeypatch
    ):
        import refgenie.managers.sources.client as client_module

        def explode(*args, **kwargs):
            raise AssertionError("a rich Progress must not be built when a sink is active")

        monkeypatch.setattr(client_module, "Progress", explode)

        events = []
        out = tmp_path / "out.bin"
        with progress.use_sink(events.append):
            download_client.download_with_progress(
                "download_op", out, url_format_params={"id": "x"}, name="thing"
            )

        assert out.read_bytes() == b"x" * 5000
        # Opens with a coarse phase announcement, then byte counts; every event
        # names its phase, which the manager turns into the UI's step counter.
        assert events[0].type == "stage"
        assert all(e.extra.get("phase") == "download" for e in events)
        byte_events = [e for e in events if e.type == "progress"]
        assert byte_events, "expected byte progress events"
        assert all(e.unit == "bytes" for e in byte_events)
        currents = [e.current for e in byte_events]
        assert currents == sorted(currents)
        assert byte_events[-1].current == 5000
        assert byte_events[-1].total == 5000

    def test_no_sink_still_builds_the_rich_display(self, download_client, tmp_path):
        """The CLI path is untouched (its behavior is covered in full by
        ``TestDownloadWithProgress`` in tests/managers/test_server_client.py)."""
        import refgenie.managers.sources.client as client_module

        built = []
        real = client_module.Progress

        class Spy(real):
            def __init__(self, *args, **kwargs):
                built.append(True)
                super().__init__(*args, **kwargs)

        client_module.Progress = Spy
        try:
            out = tmp_path / "out.bin"
            download_client.download_with_progress(
                "download_op", out, url_format_params={"id": "x"}, name="thing"
            )
        finally:
            client_module.Progress = real
        assert built == [True]
