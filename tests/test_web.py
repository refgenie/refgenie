"""The refgenie web UI: the whole dash surface in one module.

Consolidates the thirteen former dash/web-UI test files. Every contract the
web_ui_contracts note pins in its "Test ownership" table is asserted here at
least once:

* the library progress hook (``refgenie.progress``);
* ``JobManager`` lifecycle, cancel, coalescing and ring truncation;
* pull/build error mapping to the shared error vocabulary, as job records;
* the jobs HTTP API, including the one multiplexed SSE stream, replay and
  ``Last-Event-ID``;
* the local actions API ``/v1/actions/*``: HTTP envelope, 202 submission, the
  action-header / CORS / host-guard security stack, server-mode absence, and
  the synchronous curation endpoints;
* the localhost bridge: ``/ping``, bridge modes, the LNA preflight, and the
  cross-origin action policy;
* ``create_app(mode=...)`` construction, mode isolation, route/operationId
  uniqueness and SPA-route/API non-collision;
* SPA serving: index, hashed assets, cache headers, deep-link fallback,
  missing-bundle 503, the ``/service-info`` bootstrap key set;
* the OpenAPI <-> TypeScript wire-type drift guard;
* building from a worker thread; and an end-to-end local pull/build smoke.

Tier: everything defaults to ``unit`` (this file lives in ``tests/`` root).
Classes that build real genome folders/archives on disk carry an explicit
class-level ``pytest.mark.component`` -- there is deliberately NO module-level
component mark, which would wrongly promote the whole file.

Apps are always built through the real ``create_app`` factory (via the
``make_*_app``/``make_*_client`` helpers), never a hand-assembled ``FastAPI()``;
the one exception is the jobs router, which is included on a bare app on purpose
to prove it carries no dependency on the app factory.
"""

import inspect
import logging
import re
import threading
import time
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from tests.helpers import (
    ACTION_HEADERS,
    ASSET,
    BUILD_PARAMS,
    GENOME,
    GROUP,
    PHASE_VOCABULARY,
    PULL_PARAMS,
    act,
    build_params,
    instant_runner,
    job_result_ok,
    jobs_app,
    jobs_client,
    make_gated_runner,
    make_local_app,
    make_local_client,
    make_server_app,
    make_server_client,
    make_server_rgc,
    make_slow_runner,
    open_gate,
    preflight,
    pull_params,
    read_sse,
    requires_dash,
    requires_server,
    requires_web_assets,
    stub_rgc,
    submit,
    submit_and_wait,
    until_done,
    wait_for,
    web_stub_rgc,
)

requires_server()
requires_dash()

from fastapi import FastAPI  # noqa: E402  (must follow the extras guard)
from fastapi.testclient import TestClient  # noqa: E402

from refgenie import progress  # noqa: E402
from refgenie.cli.commands import asset as cli_asset  # noqa: E402
from refgenie.cli.commands import helpers as cli_helpers  # noqa: E402
from refgenie.cli.commands import listing as cli_listing  # noqa: E402
from refgenie.cli.dispatch import get_dispatch  # noqa: E402
from refgenie.exceptions import (  # noqa: E402
    AssetExistsError,
    MissingAliasError,
    MissingBuildInputError,
    MissingRecipeError,
    NoArchiveError,
    PullFailedError,
    PullSkipped,
)
from refgenie.logger import logger as refgenie_logger  # noqa: E402
from refgenie.progress import ProgressAborted  # noqa: E402
from refgenie.server.const import (  # noqa: E402
    APP_MODE_LOCAL,
    APP_MODE_SERVER,
    SPA_CLIENT_ROUTES,
)
from refgenie.server.errors import ErrorCode, install_error_handlers  # noqa: E402
from refgenie.server.jobs import jobs_router  # noqa: E402
from refgenie.server.jobs.manager import CancelOutcome, JobManager  # noqa: E402
from refgenie.server.jobs.schemas import (  # noqa: E402
    JobKind,
    JobStatus,
    PullJobParams,
)
from refgenie.server.main import CAPABILITY_KEYS, _capabilities, create_app  # noqa: E402
from refgenie.server.routers import shared as shared_router  # noqa: E402
from refgenie.server.routers import version4 as version4_router  # noqa: E402
from tests.test_server import (  # noqa: E402
    _assert_no_duplicate_routes,
    _assert_unique_operation_ids,
    _route_keys,
)

# ===========================================================================
# Shared job helpers and fixtures
# ===========================================================================

#: The result every fake runner returns, as the typed model the manager stores.
OK = job_result_ok()

#: The bridge's default allowlisted origin, and one that is never allowed.
PUBLIC_ORIGIN = "https://refgenie.org"
EVIL_ORIGIN = "https://evil.example"


@pytest.fixture
def gate():
    event = threading.Event()
    yield event
    event.set()


@pytest.fixture
def manager(job_manager_factory):
    """A default all-instant-runner manager, for the jobs router tests."""
    return job_manager_factory()


@pytest.fixture
def client(manager):
    with jobs_client(manager) as tc:
        yield tc


# ===========================================================================
# The library progress hook: refgenie.progress
# ===========================================================================


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
            "progress", "downloading", 5, 10, "bytes",
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
        ``TestDownloadWithProgress`` in tests/test_server_client.py)."""
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


# ===========================================================================
# JobManager lifecycle, against fake runners
# ===========================================================================


class TestJobManagerLifecycle:
    def test_submit_runs_and_succeeds(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: instant_runner})
        ref = manager.submit_pull(pull_params())
        assert ref.status == JobStatus.QUEUED
        assert ref.duplicate is False
        assert ref.links.events == "/v1/jobs/events"

        record = manager.wait(ref.job_id, timeout=5)
        assert record.status == JobStatus.SUCCEEDED
        assert record.result == OK
        assert record.started_at <= record.finished_at
        assert record.queue_position is None
        assert record.cancellable is False

    def test_every_transition_appends_a_status_event(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: instant_runner})
        ref = manager.submit_pull(pull_params())
        manager.wait(ref.job_id, timeout=5)

        events, cursor, truncated = manager.events_since(0)
        assert truncated is False
        assert [e.seq for e in events] == sorted(e.seq for e in events)
        assert cursor == events[-1].seq
        statuses = [e.status for e in events if e.type == "status"]
        assert statuses == [JobStatus.QUEUED, JobStatus.RUNNING, JobStatus.SUCCEEDED]
        assert [e.type for e in events][-1] == "done"

    @pytest.mark.parametrize(
        "exc,status,code,message",
        [
            (AssetExistsError("x"), JobStatus.FAILED, ErrorCode.ASSET_EXISTS,
             "Asset already exists. Use force to overwrite."),
            (RuntimeError("boom"), JobStatus.FAILED, ErrorCode.INTERNAL_ERROR, "boom"),
            (PullSkipped("nope"), JobStatus.CANCELLED, ErrorCode.PULL_SKIPPED, None),
        ],
    )
    def test_runner_exception_is_classified(self, job_manager_factory, exc, status, code, message):
        manager = job_manager_factory(
            runners={JobKind.PULL: make_slow_runner(open_gate(), steps=0, fail=exc)}
        )
        record = manager.wait(manager.submit_pull(pull_params()).job_id, timeout=5)
        assert record.status == status
        assert record.error.code == code
        if message is not None:
            assert record.error.message == message
        if isinstance(exc, AssetExistsError):
            assert "AssetExistsError" in record.error.detail

    def test_wait_times_out_rather_than_hanging(self, job_manager_factory):
        gate = threading.Event()
        manager = job_manager_factory(runners={JobKind.PULL: make_slow_runner(gate)})
        ref = manager.submit_pull(pull_params())
        with pytest.raises(TimeoutError):
            manager.wait(ref.job_id, timeout=0.3)
        gate.set()


class TestJobProgressAndEvents:
    def test_progress_is_denormalized_onto_the_record(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: make_slow_runner(open_gate())})
        record = manager.wait(manager.submit_pull(pull_params()).job_id, timeout=5)
        assert record.progress.phase == "download"
        assert record.progress.percent == 100.0

    def test_seq_is_global_across_jobs_and_events_carry_job_id(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: make_slow_runner(open_gate())})
        first = manager.submit_pull(pull_params(asset_group_name="a"))
        second = manager.submit_pull(pull_params(asset_group_name="b"))
        manager.wait(first.job_id, timeout=5)
        manager.wait(second.job_id, timeout=5)

        events, _, _ = manager.events_since(0)
        seqs = [e.seq for e in events]
        assert seqs == sorted(seqs)
        assert len(set(seqs)) == len(seqs)
        assert {first.job_id, second.job_id} == {e.job_id for e in events}

    def test_events_since_never_duplicates_or_skips(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: make_slow_runner(open_gate())})
        manager.wait(manager.submit_pull(pull_params()).job_id, timeout=5)

        collected = []
        cursor = 0
        for _ in range(5):
            events, cursor, _ = manager.events_since(cursor)
            collected.extend(events)
        all_events, _, _ = manager.events_since(0)
        assert [e.seq for e in collected] == [e.seq for e in all_events]

    def test_ring_truncation_is_reported(self, job_manager_factory):
        manager = job_manager_factory(
            runners={JobKind.PULL: make_slow_runner(open_gate(), steps=10)}, event_buffer=5
        )
        manager.wait(manager.submit_pull(pull_params()).job_id, timeout=5)
        events, _, truncated = manager.events_since(0)
        assert truncated is True
        assert len(events) == 5
        assert events[0].seq > 1


class TestJobQueueing:
    def test_queue_position_is_per_executor(self, job_manager_factory):
        """Builds and pulls are two queues. A build waiting behind a build says
        so; a pull that started in the second pull slot says nothing."""
        gate = threading.Event()
        runner = make_slow_runner(gate, steps=1)
        manager = job_manager_factory(runners={JobKind.BUILD: runner, JobKind.PULL: runner})

        first_build = manager.submit_build(build_params(asset_group_name="one"))
        wait_for(lambda: manager.record(first_build.job_id).status == JobStatus.RUNNING)
        second_build = manager.submit_build(build_params(asset_group_name="two"))
        pull = manager.submit_pull(pull_params())
        wait_for(lambda: manager.record(pull.job_id).status == JobStatus.RUNNING)

        assert manager.record(second_build.job_id).queue_position == 0
        assert manager.record(pull.job_id).queue_position is None
        gate.set()

    def test_pulls_run_concurrently(self, job_manager_factory):
        gate = threading.Event()
        manager = job_manager_factory(
            runners={JobKind.PULL: make_slow_runner(gate, steps=1)}, pull_workers=2
        )
        a = manager.submit_pull(pull_params(asset_group_name="a"))
        b = manager.submit_pull(pull_params(asset_group_name="b"))
        wait_for(
            lambda: manager.record(a.job_id).status == JobStatus.RUNNING
            and manager.record(b.job_id).status == JobStatus.RUNNING
        )
        gate.set()


class TestJobCancellation:
    def test_queued_job_is_cancelled_without_running(self, job_manager_factory):
        gate = threading.Event()
        invoked = []

        def counting_runner(ctx):
            invoked.append(ctx.job_id)
            gate.wait(timeout=5)
            return None

        manager = job_manager_factory(runners={JobKind.BUILD: counting_runner})
        first = manager.submit_build(build_params(asset_group_name="one"))
        wait_for(lambda: manager.record(first.job_id).status == JobStatus.RUNNING)
        queued = manager.submit_build(build_params(asset_group_name="two"))

        assert manager.request_cancel(queued.job_id) == CancelOutcome.ACCEPTED
        gate.set()
        record = manager.wait(queued.job_id, timeout=5)
        assert record.status == JobStatus.CANCELLED
        assert record.error.code == ErrorCode.CANCELLED
        assert queued.job_id not in invoked

    def test_running_pull_is_cancelled_at_its_next_report(self, job_manager_factory):
        gate = threading.Event()
        manager = job_manager_factory(runners={JobKind.PULL: make_slow_runner(gate, steps=50)})
        ref = manager.submit_pull(pull_params())
        wait_for(lambda: manager.record(ref.job_id).status == JobStatus.RUNNING)

        assert manager.record(ref.job_id).cancellable is True
        assert manager.request_cancel(ref.job_id) == CancelOutcome.ACCEPTED
        gate.set()
        record = manager.wait(ref.job_id, timeout=5)
        assert record.status == JobStatus.CANCELLED
        assert record.error.code == ErrorCode.CANCELLED

    def test_running_build_is_not_cancellable(self, job_manager_factory):
        gate = threading.Event()
        manager = job_manager_factory(runners={JobKind.BUILD: make_slow_runner(gate, steps=1)})
        ref = manager.submit_build(build_params())
        wait_for(lambda: manager.record(ref.job_id).status == JobStatus.RUNNING)

        assert manager.record(ref.job_id).cancellable is False
        assert manager.request_cancel(ref.job_id) == CancelOutcome.NOT_CANCELLABLE
        gate.set()

    def test_terminal_job_reports_already_done(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: instant_runner})
        ref = manager.submit_pull(pull_params())
        manager.wait(ref.job_id, timeout=5)
        assert manager.request_cancel(ref.job_id) == CancelOutcome.ALREADY_DONE


class TestJobCoalescing:
    def test_identical_in_flight_params_return_the_same_job(self, job_manager_factory):
        gate = threading.Event()
        manager = job_manager_factory(runners={JobKind.PULL: make_slow_runner(gate, steps=1)})
        first = manager.submit_pull(pull_params())
        second = manager.submit_pull(pull_params())

        assert second.job_id == first.job_id
        assert second.duplicate is True
        gate.set()
        manager.wait(first.job_id, timeout=5)
        assert len(manager.list()) == 1

    def test_a_finished_job_does_not_block_a_new_identical_one(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: instant_runner})
        first = manager.submit_pull(pull_params())
        manager.wait(first.job_id, timeout=5)
        second = manager.submit_pull(pull_params())
        assert second.job_id != first.job_id
        assert second.duplicate is False

    def test_different_params_are_different_jobs(self, job_manager_factory):
        gate = threading.Event()
        manager = job_manager_factory(runners={JobKind.PULL: make_slow_runner(gate, steps=1)})
        a = manager.submit_pull(pull_params(asset_group_name="a"))
        b = manager.submit_pull(pull_params(asset_group_name="b"))
        assert a.job_id != b.job_id
        gate.set()


class TestJobHistoryAndShutdown:
    def test_history_evicts_only_finished_jobs_oldest_first(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: instant_runner}, history_limit=3)
        ids = []
        for i in range(5):
            ref = manager.submit_pull(pull_params(asset_group_name=f"g{i}"))
            manager.wait(ref.job_id, timeout=5)
            ids.append(ref.job_id)
        remaining = {record.id for record in manager.list()}
        assert remaining == set(ids[-3:])

    def test_running_jobs_are_never_evicted(self, job_manager_factory):
        gate = threading.Event()
        manager = job_manager_factory(
            runners={JobKind.PULL: instant_runner, JobKind.BUILD: make_slow_runner(gate, steps=1)},
            history_limit=2,
        )
        running = manager.submit_build(build_params())
        wait_for(lambda: manager.record(running.job_id).status == JobStatus.RUNNING)
        for i in range(5):
            ref = manager.submit_pull(pull_params(asset_group_name=f"g{i}"))
            manager.wait(ref.job_id, timeout=5)
        assert manager.record(running.job_id).status == JobStatus.RUNNING
        gate.set()

    def test_forget_drops_a_terminal_job_and_refuses_a_running_one(self, job_manager_factory):
        gate = threading.Event()
        manager = job_manager_factory(
            runners={JobKind.PULL: instant_runner, JobKind.BUILD: make_slow_runner(gate, steps=1)}
        )
        done = manager.submit_pull(pull_params())
        manager.wait(done.job_id, timeout=5)
        manager.forget(done.job_id)
        with pytest.raises(KeyError):
            manager.record(done.job_id)

        running = manager.submit_build(build_params())
        wait_for(lambda: manager.record(running.job_id).status == JobStatus.RUNNING)
        with pytest.raises(ValueError):
            manager.forget(running.job_id)
        gate.set()

    def test_shutdown_returns_promptly_with_a_job_running(self, job_manager_factory, caplog):
        gate = threading.Event()
        manager = job_manager_factory(runners={JobKind.BUILD: make_slow_runner(gate, steps=1)})
        ref = manager.submit_build(build_params())
        wait_for(lambda: manager.record(ref.job_id).status == JobStatus.RUNNING)

        with caplog.at_level(logging.WARNING, logger="refgenie"):
            started = time.monotonic()
            manager.shutdown(wait=False)
            elapsed = time.monotonic() - started
        assert elapsed < 1.0
        assert "abandoned at process exit" in caplog.text
        gate.set()

    def test_listing_filters_by_kind_and_status(self, job_manager_factory):
        manager = job_manager_factory(
            runners={JobKind.PULL: instant_runner, JobKind.BUILD: instant_runner}
        )
        pull = manager.submit_pull(pull_params())
        build = manager.submit_build(build_params())
        manager.wait(pull.job_id, timeout=5)
        manager.wait(build.job_id, timeout=5)

        assert [r.id for r in manager.list(kind=JobKind.PULL)] == [pull.job_id]
        assert len(manager.list(status=JobStatus.SUCCEEDED)) == 2
        assert manager.list(limit=1)[0].id == build.job_id  # newest first


class TestJobLogCapture:
    def test_log_lines_reach_the_right_job(self, job_manager_factory):
        gate = threading.Event()

        def logging_runner(ctx):
            refgenie_logger.info(f"hello from {ctx.params.asset_group_name}")
            gate.wait(timeout=5)
            return None

        manager = job_manager_factory(runners={JobKind.PULL: logging_runner}, pull_workers=2)
        a = manager.submit_pull(pull_params(asset_group_name="aaa"))
        b = manager.submit_pull(pull_params(asset_group_name="bbb"))
        gate.set()
        manager.wait(a.job_id, timeout=5)
        manager.wait(b.job_id, timeout=5)

        events, _, _ = manager.events_since(0)
        logs = {e.job_id: [] for e in events}
        for event in events:
            if event.type == "log":
                assert event.source == "refgenie"
                logs[event.job_id].append(event.line)
        assert "hello from aaa" in logs[a.job_id]
        assert "hello from bbb" in logs[b.job_id]
        assert not any("bbb" in m for m in logs[a.job_id])

    def test_a_broken_log_route_cannot_kill_the_job(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: _logging_instant_runner})
        manager._route_log = lambda *args: 1 / 0  # sabotage the route the handler calls
        manager._log_handler._route = manager._route_log
        record = manager.wait(manager.submit_pull(pull_params()).job_id, timeout=5)
        assert record.status == JobStatus.SUCCEEDED

    def test_logs_from_outside_a_job_are_ignored(self, job_manager_factory):
        manager = job_manager_factory(runners={JobKind.PULL: instant_runner})
        ref = manager.submit_pull(pull_params())
        manager.wait(ref.job_id, timeout=5)
        before, _, _ = manager.events_since(0)
        refgenie_logger.info("a message from the main thread")
        after, _, _ = manager.events_since(0)
        assert len(after) == len(before)


def _logging_instant_runner(ctx):
    refgenie_logger.info("working")
    return OK


# ===========================================================================
# Pull and build error mapping, as job records (ex-TestDashPullEndpointErrors)
# ===========================================================================

RUNNER_PULL = pull_params(
    server_url="http://mock-server", genome_name="test_genome", force=True
)
RUNNER_BUILD = build_params(genome_name="test_genome")


def fake_asset(digest="a" * 64, name="default", registry_path="dig/fasta:default"):
    """A stand-in for a returned Asset row, including the relationship chain
    ``registry_path`` walks."""
    asset = MagicMock()
    asset.digest = digest
    asset.name = name
    asset.registry_path = registry_path
    asset.asset_group.name = "fasta"
    asset.asset_group.genome.digest = "g" * 64
    return asset


@pytest.fixture
def run_job():
    """Run one real runner against a mocked Refgenie and return the record."""
    managers = []

    def run(kind, params, rgc):
        manager = JobManager(rgc, history_limit=10)
        managers.append(manager)
        submit_fn = {"pull": manager.submit_pull, "build": manager.submit_build}[kind]
        ref = submit_fn(params)
        return manager.wait(ref.job_id, timeout=10)

    yield run
    for manager in managers:
        manager.shutdown(wait=False)


class TestPullErrors:
    """The mapping table, one row per case."""

    @pytest.mark.parametrize(
        "exc,expected_status,expected_code",
        [
            (AssetExistsError("test/fasta:default already exists"), JobStatus.FAILED, ErrorCode.ASSET_EXISTS),
            (NoArchiveError("No archive found for asset digest abc123"), JobStatus.FAILED, ErrorCode.NO_ARCHIVE),
            (PullFailedError("Server refused"), JobStatus.FAILED, ErrorCode.PULL_FAILED),
            (PullSkipped("Skipping pull of x"), JobStatus.CANCELLED, ErrorCode.PULL_SKIPPED),
            (RuntimeError("something unexpected"), JobStatus.FAILED, ErrorCode.INTERNAL_ERROR),
        ],
    )
    def test_exception_maps_to_status_and_code(self, run_job, exc, expected_status, expected_code):
        rgc = MagicMock()
        rgc.pull.side_effect = exc
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.status == expected_status
        assert record.error.code == expected_code
        assert record.error.detail, "a traceback belongs on the record"

    def test_asset_exists_message_is_actionable(self, run_job):
        rgc = MagicMock()
        rgc.pull.side_effect = AssetExistsError("test/fasta:default already exists")
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.error.message == "Asset already exists. Use force to overwrite."

    @pytest.mark.parametrize(
        "exc,fragment",
        [
            (NoArchiveError("No archive found for asset digest abc123"), "no archive found"),
            (PullFailedError("No server subscriptions found"), "no server subscriptions"),
        ],
    )
    def test_the_actual_message_survives(self, run_job, exc, fragment):
        rgc = MagicMock()
        rgc.pull.side_effect = exc
        record = run_job("pull", RUNNER_PULL, rgc)
        assert fragment in record.error.message.lower()

    def test_returning_none_is_no_subscriptions(self, run_job):
        """``pull`` returns None when nothing is subscribed and the subscribe
        prompt was declined -- which, from a job, it always is."""
        rgc = MagicMock()
        rgc.pull.return_value = None
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.status == JobStatus.FAILED
        assert record.error.code == ErrorCode.NO_SUBSCRIPTIONS
        assert "subscribe" in record.error.message.lower()


class TestPullSuccess:
    def test_result_carries_what_the_ui_needs_to_refresh(self, run_job):
        rgc = MagicMock()
        rgc.pull.return_value = fake_asset()
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.status == JobStatus.SUCCEEDED
        assert record.result.asset_digest == "a" * 64
        assert record.result.registry_path == "dig/fasta:default"
        assert record.result.asset_group_name == "fasta"
        assert record.result.genome_digest == "g" * 64
        assert record.result.staged is False

    def test_the_facade_is_called_non_interactively(self, run_job):
        """``force_large=True`` because a browser cannot answer a size prompt,
        and an explicit ``confirm`` because ``resolve_confirmer(None)`` falls
        through to a real ``rich`` prompt once the CLI has enabled prompts."""
        rgc = MagicMock()
        rgc.pull.return_value = fake_asset()
        run_job("pull", RUNNER_PULL, rgc)

        kwargs = rgc.pull.call_args.kwargs
        assert kwargs["force_large"] is True
        assert kwargs["confirm"] is not None
        assert kwargs["confirm"]("replace everything?") is False
        assert kwargs["force_server_urls"] == ["http://mock-server"]
        assert kwargs["alias_name"] == "test_genome"
        assert kwargs["asset_group_name"] == "fasta"
        assert kwargs["force"] is True

    def test_a_detached_asset_does_not_fail_the_job(self, run_job):
        """A successful pull must not be reported as a failure because its
        receipt could not be formatted."""
        asset = MagicMock()
        asset.digest = "a" * 64
        asset.name = "default"
        type(asset).registry_path = property(
            lambda self: (_ for _ in ()).throw(RuntimeError("detached"))
        )
        rgc = MagicMock()
        rgc.pull.return_value = asset
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.status == JobStatus.SUCCEEDED
        assert record.result.asset_digest == "a" * 64
        assert record.result.registry_path is None


class TestBuildErrors:
    def test_success(self, run_job, tmp_path):
        rgc = MagicMock()
        rgc.genome_folder = tmp_path
        rgc.build_asset.return_value = fake_asset()
        record = run_job("build", RUNNER_BUILD, rgc)
        assert record.status == JobStatus.SUCCEEDED
        assert record.result.registry_path == "dig/fasta:default"
        assert record.result.asset_digest == "a" * 64

    def test_none_return_is_build_failed(self, run_job, tmp_path):
        """``build_asset`` returning None means the pipeline failed. The library
        answer is a printed message and exit 0; a job card needs a code."""
        rgc = MagicMock()
        rgc.genome_folder = tmp_path
        rgc.build_asset.return_value = None
        record = run_job("build", RUNNER_BUILD, rgc)
        assert record.status == JobStatus.FAILED
        assert record.error.code == ErrorCode.BUILD_FAILED
        assert "pipeline log" in record.error.message

    def test_the_pypiper_log_is_recorded_on_the_job(self, run_job, tmp_path):
        """The tailer latches onto the build's log file, so a failed build can
        be linked to its output."""
        log_dir = tmp_path / "builds" / "test_genome" / "fasta"
        log_dir.mkdir(parents=True)
        log_file = log_dir / "refgenie_test_genome_fasta_default_log.md"

        rgc = MagicMock()
        rgc.genome_folder = tmp_path

        def slow_build(**kwargs):
            log_file.write_text("### Pipeline started\ncommand: samtools faidx\n")
            time.sleep(1.2)
            return None

        rgc.build_asset.side_effect = slow_build
        record = run_job("build", RUNNER_BUILD, rgc)
        assert record.status == JobStatus.FAILED
        assert record.log_file == str(log_file)

    @pytest.mark.parametrize(
        "exc,expected_code",
        [
            (MissingRecipeError("fasta"), ErrorCode.RECIPE_NOT_FOUND),
            (MissingAliasError("rCRSd"), ErrorCode.ALIAS_NOT_FOUND),
            (MissingBuildInputError("needs a fasta"), ErrorCode.MISSING_BUILD_INPUT),
            (RuntimeError("boom"), ErrorCode.INTERNAL_ERROR),
        ],
    )
    def test_build_exceptions_use_the_shared_vocabulary(self, run_job, tmp_path, exc, expected_code):
        rgc = MagicMock()
        rgc.genome_folder = tmp_path
        rgc.build_asset.side_effect = exc
        record = run_job("build", RUNNER_BUILD, rgc)
        assert record.status == JobStatus.FAILED
        assert record.error.code == expected_code


# ===========================================================================
# The jobs HTTP API, including the multiplexed SSE stream
#
# The app here is a bare FastAPI() with the real router included -- deliberate:
# the contract is that the jobs router carries no dependency on the app factory.
# Its mode isolation *inside* create_app is asserted in the app-modes tests.
# ===========================================================================


class TestJobsSubmission:
    def test_post_returns_201_a_job_ref_and_a_location(self, client, manager):
        response = client.post("/v1/jobs", json={"kind": "pull", "params": PULL_PARAMS})
        assert response.status_code == 201
        body = response.json()
        assert body["kind"] == "pull"
        assert body["status"] in ("queued", "running", "succeeded")
        assert body["duplicate"] is False
        assert body["links"]["events"] == "/v1/jobs/events"
        assert response.headers["Location"] == body["links"]["self"]

        record = manager.wait(body["job_id"], timeout=5)
        assert record.params["asset_group_name"] == "fasta"

    def test_pull_without_a_genome_is_422(self, client):
        response = client.post(
            "/v1/jobs", json={"kind": "pull", "params": {"asset_group_name": "fasta"}}
        )
        assert response.status_code == 422

    def test_unknown_kind_is_422(self, client):
        assert client.post("/v1/jobs", json={"kind": "nope", "params": {}}).status_code == 422

    def test_unknown_param_is_rejected(self, client):
        response = client.post(
            "/v1/jobs", json={"kind": "pull", "params": {**PULL_PARAMS, "typo": 1}}
        )
        assert response.status_code == 422

    def test_staging_without_a_stage_folder_is_400_and_creates_no_job(self, client, manager):
        response = client.post(
            "/v1/jobs", json={"kind": "build", "params": {**BUILD_PARAMS, "stage": True}}
        )
        assert response.status_code == 400
        assert "genome_stage_folder" in response.json()["detail"]
        assert manager.list() == []

    def test_staging_with_a_stage_folder_is_accepted(self, job_manager_factory, tmp_path):
        manager = job_manager_factory()
        manager.refgenie.genome_stage_folder = tmp_path
        with jobs_client(manager) as tc:
            response = tc.post(
                "/v1/jobs", json={"kind": "build", "params": {**BUILD_PARAMS, "stage": True}}
            )
        assert response.status_code == 201


class TestJobsReads:
    def test_get_job_and_404(self, client):
        job_id = submit(client)["job_id"]
        assert client.get(f"/v1/jobs/{job_id}").status_code == 200
        assert client.get("/v1/jobs/nosuchjob").status_code == 404

    def test_list_filters(self, client, manager):
        pull = submit(client)
        build = submit(client, "build")
        manager.wait(pull["job_id"], timeout=5)
        manager.wait(build["job_id"], timeout=5)

        jobs = client.get("/v1/jobs", params={"kind": "pull"}).json()["items"]
        assert [j["id"] for j in jobs] == [pull["job_id"]]
        assert len(client.get("/v1/jobs", params={"status": "succeeded"}).json()["items"]) == 2
        assert client.get("/v1/jobs", params={"status": "running"}).json()["items"] == []

    def test_polling_fallback_returns_this_jobs_events_and_the_record(self, client, manager):
        first = submit(client)
        second = submit(client, "build")
        manager.wait(first["job_id"], timeout=5)
        manager.wait(second["job_id"], timeout=5)

        body = client.get(f"/v1/jobs/{first['job_id']}/events").json()
        assert {e["job_id"] for e in body["events"]} == {first["job_id"]}
        assert body["job"]["id"] == first["job_id"]
        assert body["next"] >= body["events"][-1]["seq"]
        assert body["truncated"] is False

        resumed = client.get(
            f"/v1/jobs/{first['job_id']}/events", params={"since": body["next"]}
        ).json()
        assert resumed["events"] == []

    def test_polling_fallback_404s_on_an_unknown_job(self, client):
        assert client.get("/v1/jobs/nope/events").status_code == 404

    def test_log_endpoint(self, client, manager, tmp_path):
        job_id = submit_and_wait(client, manager)

        empty = client.get(f"/v1/jobs/{job_id}/log").json()
        assert empty == {"lines": [], "next_offset": 0, "truncated": False}

        log = tmp_path / "pipeline_log.md"
        log.write_text("one\ntwo\nthree\n")
        manager.get(job_id).log_file = str(log)

        first_page = client.get(f"/v1/jobs/{job_id}/log", params={"limit": 2}).json()
        assert first_page == {"lines": ["one", "two"], "next_offset": 2, "truncated": True}
        second_page = client.get(f"/v1/jobs/{job_id}/log", params={"offset": 2}).json()
        assert second_page == {"lines": ["three"], "next_offset": 3, "truncated": False}


class TestJobsCancelAndForget:
    def test_cancel_a_running_pull(self, job_manager_factory, gate):
        manager = job_manager_factory(runners={JobKind.PULL: make_gated_runner(gate)})
        with jobs_client(manager) as tc:
            job_id = submit(tc)["job_id"]
            wait_for(lambda: manager.record(job_id).status == JobStatus.RUNNING)
            response = tc.post(f"/v1/jobs/{job_id}/cancel")
        assert response.status_code == 202
        assert response.json()["id"] == job_id

    def test_cancel_a_running_build_is_409_not_cancellable(self, job_manager_factory, gate):
        manager = job_manager_factory(runners={JobKind.BUILD: make_gated_runner(gate)})
        with jobs_client(manager) as tc:
            job_id = submit(tc, "build")["job_id"]
            wait_for(lambda: manager.record(job_id).status == JobStatus.RUNNING)
            response = tc.post(f"/v1/jobs/{job_id}/cancel")
        assert response.status_code == 409
        assert response.json()["detail"]["code"] == "not_cancellable"

    def test_cancel_a_finished_job_is_409_already_done(self, client, manager):
        job_id = submit_and_wait(client, manager)
        response = client.post(f"/v1/jobs/{job_id}/cancel")
        assert response.status_code == 409
        assert response.json()["detail"]["code"] == "already_done"

    def test_cancel_an_unknown_job_is_404(self, client):
        assert client.post("/v1/jobs/nope/cancel").status_code == 404

    def test_delete_forgets_a_terminal_job_and_refuses_a_running_one(self, job_manager_factory, gate):
        manager = job_manager_factory(
            runners={JobKind.PULL: instant_runner, JobKind.BUILD: make_gated_runner(gate)}
        )
        with jobs_client(manager) as tc:
            done = submit(tc)
            manager.wait(done["job_id"], timeout=5)
            assert tc.delete(f"/v1/jobs/{done['job_id']}").status_code == 204
            assert tc.get(f"/v1/jobs/{done['job_id']}").status_code == 404

            running = submit(tc, "build")
            wait_for(lambda: manager.record(running["job_id"]).status == JobStatus.RUNNING)
            assert tc.delete(f"/v1/jobs/{running['job_id']}").status_code == 409


class TestJobsServerMode:
    """No manager on app.state -- which is what server mode looks like."""

    @pytest.fixture
    def bare_client(self):
        app = FastAPI()
        app.include_router(jobs_router, prefix="/v1")
        with TestClient(app, headers=ACTION_HEADERS) as tc:
            yield tc

    @pytest.mark.parametrize(
        "method,path",
        [
            ("get", "/v1/jobs"),
            ("get", "/v1/jobs/abc"),
            ("get", "/v1/jobs/abc/events"),
            ("get", "/v1/jobs/abc/log"),
            ("get", "/v1/jobs/events"),
            ("post", "/v1/jobs/abc/cancel"),
            ("delete", "/v1/jobs/abc"),
        ],
    )
    def test_every_endpoint_is_503(self, bare_client, method, path):
        assert getattr(bare_client, method)(path).status_code == 503


class TestJobsWriteGuards:
    """The three write routes carry the same ``X-Refgenie-Action`` guard as the
    actions router -- attached per-route, because the reads and the SSE stream
    (which ``EventSource`` cannot send headers to) must stay open. The
    full-stack 403 is owned by the actions security tests; here the real error
    handlers are installed so the envelope shape is asserted too."""

    WRITE_ROUTES = [
        ("post", "/v1/jobs", {"kind": "pull", "params": PULL_PARAMS}),
        ("post", "/v1/jobs/abc/cancel", None),
        ("delete", "/v1/jobs/abc", None),
    ]

    @pytest.fixture
    def guarded_client(self, manager):
        app = jobs_app(manager)
        install_error_handlers(app)
        with TestClient(app) as tc:  # no default action header, on purpose
            yield tc

    @pytest.mark.parametrize("method,path,body", WRITE_ROUTES)
    def test_writes_without_the_header_are_403_with_envelope(
        self, guarded_client, manager, method, path, body
    ):
        response = getattr(guarded_client, method)(
            path, **({"json": body} if body is not None else {})
        )
        assert response.status_code == 403, f"{method} {path}"
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "missing_action_header"
        assert manager.list() == []  # nothing was submitted

    @pytest.mark.parametrize("method,path,body", WRITE_ROUTES)
    def test_writes_with_the_header_reach_the_handler(self, guarded_client, method, path, body):
        response = getattr(guarded_client, method)(
            path, headers=ACTION_HEADERS, **({"json": body} if body is not None else {})
        )
        assert response.status_code != 403, f"{method} {path}: {response.text}"

    def test_reads_and_the_stream_require_no_header(self, guarded_client, manager):
        job_id = guarded_client.post(
            "/v1/jobs", json={"kind": "pull", "params": PULL_PARAMS}, headers=ACTION_HEADERS
        ).json()["job_id"]
        manager.wait(job_id, timeout=5)
        for path in ("/v1/jobs", f"/v1/jobs/{job_id}", f"/v1/jobs/{job_id}/events",
                     f"/v1/jobs/{job_id}/log"):
            assert guarded_client.get(path).status_code == 200, path


class TestJobsSSE:
    """The one multiplexed stream. These drive the ASGI app directly instead of
    ``TestClient``, which buffers the whole body and can never read a stream
    that stays open -- the defining property of the endpoint under test."""

    def test_replays_after_completion_and_stays_open(self, job_manager_factory):
        """A stream opened AFTER a job finished still replays it and still
        delivers ``done`` -- that removes the submit-then-connect race."""
        manager = job_manager_factory()
        first = manager.submit_pull(pull_params())
        manager.wait(first.job_id, timeout=5)

        status, headers, frames = read_sse(manager, until=until_done)

        assert status == 200
        assert headers["content-type"].startswith("text/event-stream")
        assert headers["cache-control"] == "no-cache"
        assert headers["x-accel-buffering"] == "no"

        names = [f["event"] for f in frames]
        assert names[0] == "status"
        assert names[-1] == "done"
        done = frames[-1]["data"]
        assert done["job_id"] == first.job_id
        assert done["id"] == first.job_id
        assert done["status"] == "succeeded"
        assert done["result"]["asset_digest"] == "d" * 64
        assert done["links"]["self"].endswith(first.job_id)
        assert [int(f["id"]) for f in frames] == sorted(int(f["id"]) for f in frames)

    def test_stream_carries_a_later_job_on_the_same_connection(self, job_manager_factory):
        """The stream does not end at the first ``done``: it is multiplexed."""
        manager = job_manager_factory()
        first = manager.submit_pull(pull_params())
        manager.wait(first.job_id, timeout=5)
        later = {}

        async def submit_a_second_job():
            import asyncio

            await asyncio.sleep(0.3)
            later["ref"] = manager.submit_build(build_params())

        _, _, frames = read_sse(
            manager,
            until=lambda f: f["event"] == "done"
            and later.get("ref") is not None
            and f["data"]["job_id"] == later["ref"].job_id,
            during=submit_a_second_job,
        )
        done_ids = [f["data"]["job_id"] for f in frames if f["event"] == "done"]
        assert done_ids == [first.job_id, later["ref"].job_id]

    def test_last_event_id_and_since_resume_from_the_global_sequence(self, job_manager_factory):
        manager = job_manager_factory()
        job = manager.submit_pull(pull_params())
        manager.wait(job.job_id, timeout=5)

        _, _, by_header = read_sse(
            manager, until=until_done, headers=[(b"last-event-id", b"2")]
        )
        assert int(by_header[0]["id"]) == 3

        _, _, by_query = read_sse(manager, until=until_done, query="since=2")
        assert [f["id"] for f in by_query] == [f["id"] for f in by_header]

    def test_no_custom_header_is_required(self, job_manager_factory):
        """``EventSource`` cannot set headers, so the action-header dependency
        must never be attached to this route."""
        manager = job_manager_factory()
        job = manager.submit_pull(pull_params())
        manager.wait(job.job_id, timeout=5)
        status, _, _ = read_sse(manager, until=until_done, headers=[])
        assert status == 200

    def test_keepalive_is_a_named_heartbeat_event(self, job_manager_factory, monkeypatch):
        """Not a ``: keepalive`` comment -- EventSource never dispatches those
        to JS, so a client watchdog would see nothing and reconnect forever."""
        import refgenie.server.jobs.router as router_module

        monkeypatch.setattr(router_module, "HEARTBEAT_SECONDS", 0.25)
        manager = job_manager_factory()
        _, _, frames = read_sse(manager, until=lambda f: f["event"] == "heartbeat")
        assert frames[-1]["event"] == "heartbeat"
        assert "seq" in frames[-1]["data"]

    def test_truncation_is_announced(self, job_manager_factory):
        manager = job_manager_factory(event_buffer=2)
        job = manager.submit_pull(pull_params())
        manager.wait(job.job_id, timeout=5)
        _, _, frames = read_sse(manager, until=until_done)
        assert frames[0]["event"] == "truncated"
        assert frames[0]["data"]["dropped"] >= 1


class TestJobsWireContract:
    """Field-for-field agreement with ``frontend/src/services/contracts.ts``.
    A rename here is not a compatible change: the client keys a map on
    ``undefined``, merges nulls into its progress object, or -- for an event
    name it does not know -- never receives the frame at all."""

    def test_a_record_is_keyed_id_and_a_ref_is_keyed_job_id(self, client, manager):
        """``JobRef`` (a submission receipt) is keyed ``job_id``; ``Job`` (job
        state) is keyed ``id``. The client's map is keyed on ``id`` and its
        "record or patch?" test is ``typeof candidate.id === 'string'``."""
        ref = submit(client)
        assert "job_id" in ref and "id" not in ref

        record = client.get(f"/v1/jobs/{ref['job_id']}").json()
        assert record["id"] == ref["job_id"]
        assert "job_id" not in record

    def test_a_record_carries_every_field_the_console_renders(self, client, manager):
        job_id = submit_and_wait(client, manager)
        record = client.get(f"/v1/jobs/{job_id}").json()

        for field in (
            "id", "kind", "status", "label", "target", "queue_position", "progress",
            "created_at", "started_at", "finished_at", "result", "error",
            "log_lines", "cancellable",
        ):
            assert field in record, f"the client reads {field!r}"

        assert record["label"] == "pull rCRSd/fasta"
        assert record["target"] == {
            "genome_digest": None,
            "genome_name": "rCRSd",
            "asset_group_name": "fasta",
            "asset_name": None,
        }
        assert record["result"]["asset_digest"] == "d" * 64

    def test_list_uses_the_same_envelope_as_every_other_collection(self, client, manager):
        job_id = submit_and_wait(client, manager)

        body = client.get("/v1/jobs").json()
        assert set(body) == {"items", "pagination"}
        assert set(body["pagination"]) == {"offset", "limit", "total"}
        assert body["items"][0]["id"] == job_id

    def test_status_accepts_the_active_and_terminal_groupings(self, job_manager_factory, gate):
        """The console asks these two questions on every poll cycle."""
        manager = job_manager_factory(
            runners={JobKind.PULL: instant_runner, JobKind.BUILD: make_gated_runner(gate)}
        )
        with jobs_client(manager) as tc:
            done = submit(tc)
            manager.wait(done["job_id"], timeout=5)
            running = submit(tc, "build")
            wait_for(lambda: manager.record(running["job_id"]).status == JobStatus.RUNNING)

            active = tc.get("/v1/jobs", params={"status": "active"}).json()["items"]
            terminal = tc.get("/v1/jobs", params={"status": "terminal"}).json()["items"]
            assert [j["id"] for j in active] == [running["job_id"]]
            assert [j["id"] for j in terminal] == [done["job_id"]]
            assert tc.get("/v1/jobs", params={"status": "nonsense"}).status_code == 422

    def test_pagination_windows_the_list(self, client, manager):
        for index in range(3):
            ref = client.post(
                "/v1/jobs",
                json={"kind": "pull", "params": {**PULL_PARAMS, "asset_group_name": f"g{index}"}},
            ).json()
            manager.wait(ref["job_id"], timeout=5)

        body = client.get("/v1/jobs", params={"offset": 1, "limit": 1}).json()
        assert len(body["items"]) == 1
        assert body["pagination"] == {"offset": 1, "limit": 1, "total": 3}

    def test_only_the_six_known_event_names_are_emitted(self, job_manager_factory):
        """A seventh name is not an extension -- the browser drops it."""
        manager = job_manager_factory()
        job = manager.submit_pull(pull_params())
        manager.wait(job.job_id, timeout=5)
        _, _, frames = read_sse(manager, until=until_done)
        assert {f["event"] for f in frames} <= {
            "status", "progress", "log", "done", "truncated", "heartbeat"
        }

    def test_a_progress_frame_carries_the_flattened_progress_fields(self, job_manager_factory):
        """The client destructures the frame and merges the remainder, so the
        names must be the ``JobProgress`` ones and they must be top-level."""
        manager = job_manager_factory(runners={JobKind.PULL: _byte_progress_runner})
        job = manager.submit_pull(pull_params())
        manager.wait(job.job_id, timeout=5)

        _, _, frames = read_sse(manager, until=until_done)
        progress_frames = [f["data"] for f in frames if f["event"] == "progress"]
        assert progress_frames, "expected at least one progress frame"
        byte_frame = progress_frames[-1]
        assert byte_frame["phase"] == "download"
        assert byte_frame["bytes_done"] == 512
        assert byte_frame["bytes_total"] == 1024
        assert byte_frame["percent"] == 50.0
        assert "line" not in byte_frame and "status" not in byte_frame

    def test_a_log_frame_carries_line_and_source_not_message(self, job_manager_factory):
        """``message`` is not read on a log frame; the client does
        ``lines ?? [line]``, so a frame with neither is dropped."""
        manager = job_manager_factory(runners={JobKind.PULL: _logging_runner})
        job = manager.submit_pull(pull_params())
        manager.wait(job.job_id, timeout=5)

        _, _, frames = read_sse(manager, until=until_done)
        logs = [f["data"] for f in frames if f["event"] == "log"]
        assert logs, "expected a log frame"
        assert logs[0]["line"] == "a line from the runner"
        assert logs[0]["source"] == "refgenie"
        assert "message" not in logs[0]

    def test_a_queued_status_frame_reports_its_queue_position(self, job_manager_factory, gate):
        manager = job_manager_factory(runners={JobKind.BUILD: make_gated_runner(gate)})
        first = manager.submit_build(build_params())
        wait_for(lambda: manager.record(first.job_id).status == JobStatus.RUNNING)
        manager.submit_build(build_params(asset_group_name="second"))

        _, _, frames = read_sse(
            manager,
            until=lambda f: f["event"] == "status" and f["data"].get("queue_position") == 0,
        )
        assert frames[-1]["data"]["status"] == "queued"


def _byte_progress_runner(ctx):
    progress.emit(
        "progress", message="thing", current=512, total=1024, unit="bytes", phase="download"
    )
    return OK


def _logging_runner(ctx):
    ctx.log("a line from the runner")
    return OK



# ===========================================================================
# The local actions API: /v1/actions/* (refgenie/server/local/)
# ===========================================================================

#: Every state-changing route, with a minimal valid body. The security tests
#: parametrize over this list so a newly added route is covered by adding a row.
ROUTES = [
    ("POST", "/v1/actions/pull", {"asset_group": "fasta", "genome": "rCRSd"}),
    ("POST", "/v1/actions/build", {"recipe": "fasta", "genome": "rCRSd", "asset_group": "fasta"}),
    ("POST", "/v1/actions/build/preflight",
     {"recipe": "fasta", "genome": "rCRSd", "asset_group": "fasta"}),
    ("POST", "/v1/actions/genomes", {"fasta": "/tmp/genome.fa", "aliases": ["g1"]}),
    ("DELETE", "/v1/actions/assets/deadbeef", None),
    ("DELETE", "/v1/actions/genomes/rCRSd", None),
    ("POST", "/v1/actions/aliases", {"alias": "a2", "genome_digest": "d1"}),
    ("DELETE", "/v1/actions/aliases/a2", None),
    ("POST", "/v1/actions/subscriptions", {"server_urls": ["http://s.example"]}),
    ("DELETE", "/v1/actions/subscriptions", {"server_urls": ["http://s.example"]}),
    ("POST", "/v1/actions/assets/default",
     {"genome_digest": "d1", "asset_group": "fasta", "asset": "test"}),
]


@pytest.mark.parametrize(
    "env_var,value",
    [
        ("REFGENIE_LOCAL_ALLOWED_ORIGINS", '["*"]'),
        ("REFGENIE_BRIDGE_ORIGINS", "*"),
        ("REFGENIE_BRIDGE_ORIGIN_REGEX", ".*"),
    ],
)
def test_wildcard_origin_config_refuses_to_construct(monkeypatch, env_var, value):
    """``allow_origins=['*']`` (or a ``.*`` regex) would defeat the whole design
    in local mode: fail at app-construction time, whichever knob sets it."""
    monkeypatch.setenv(env_var, value)
    with pytest.raises(ValueError, match="allowlist"):
        make_local_app(web_stub_rgc())


class TestActionSecurity:
    """Header, CORS allowlist, host guard, and the server-mode absence."""

    @pytest.fixture()
    def client(self):
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as c:
            yield c

    @pytest.mark.parametrize("method,path,body", ROUTES)
    def test_missing_action_header_is_403_with_envelope(self, client, method, path, body):
        response = client.request(method, path, json=body)
        assert response.status_code == 403, f"{method} {path}"
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "missing_action_header"

    @pytest.mark.parametrize("method,path,body", ROUTES)
    def test_with_header_the_request_reaches_the_handler(self, client, method, path, body):
        response = act(client, method, path, json=body)
        assert response.status_code != 403, f"{method} {path}: {response.text}"

    def test_preflight_for_allowed_origin_allows_the_action_header(self, client):
        response = preflight(
            client, "/v1/actions/pull", PUBLIC_ORIGIN, "POST",
            **{"Access-Control-Request-Headers": "x-refgenie-action"},
        )
        assert response.status_code == 200
        assert response.headers["access-control-allow-origin"] == PUBLIC_ORIGIN
        assert "x-refgenie-action" in response.headers["access-control-allow-headers"].lower()

    def test_preflight_for_delete_is_refused(self, client):
        """Destructive verbs are not on the cross-origin surface at all."""
        response = preflight(client, "/v1/actions/assets/deadbeef", PUBLIC_ORIGIN, "DELETE")
        assert response.status_code == 400

    def test_preflight_for_unlisted_origin_gets_no_cors_grant(self, client):
        response = preflight(client, "/v1/actions/pull", EVIL_ORIGIN, "POST")
        assert "access-control-allow-origin" not in response.headers

    def test_origin_allowlist_is_configurable(self, monkeypatch):
        # Bridge off so this isolates REFGENIE_LOCAL_ALLOWED_ORIGINS: with the
        # bridge on (default "read"), refgenie.org would be granted anyway.
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        monkeypatch.setenv("REFGENIE_LOCAL_ALLOWED_ORIGINS", '["https://example.test"]')
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            granted = preflight(client, "/v1/actions/pull", "https://example.test", "POST")
            refused = preflight(client, "/v1/actions/pull", PUBLIC_ORIGIN, "POST")
        assert granted.headers["access-control-allow-origin"] == "https://example.test"
        assert "access-control-allow-origin" not in refused.headers

    def test_non_loopback_host_is_421_forbidden_host(self, client):
        """The DNS-rebinding guard: a rebound hostname arrives as a foreign Host
        header and must be refused before any handler runs."""
        response = client.get("/service-info", headers={"host": "evil.example"})
        assert response.status_code == 421
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "forbidden_host"

    @pytest.mark.parametrize("host", ["127.0.0.1", "localhost:8080", "[::1]:9000"])
    def test_loopback_hosts_pass_the_guard(self, client, host):
        assert client.get("/service-info", headers={"host": host}).status_code == 200

    def test_server_mode_has_no_actions_routes(self):
        """Security-critical: several request models accept server-local
        filesystem paths, so the router must never reach the public server."""
        app = make_server_app(stub_rgc())
        actions_paths = [p for _, p in _route_keys(app) if p.startswith("/v1/actions")]
        assert actions_paths == []


class TestPullEndpoint:
    """The submission contract only -- error mapping lives with the runners."""

    def test_pull_returns_202_jobref_and_typed_params(self):
        rgc = web_stub_rgc()
        with make_local_client(rgc) as client:
            response = act(
                client, "POST", "/v1/actions/pull",
                json={"asset_group": "fasta", "genome": "rCRSd"},
            )
            assert response.status_code == 202
            ref = response.json()
            assert ref["kind"] == "pull"
            assert ref["status"] == "queued"
            assert ref["duplicate"] is False
            assert ref["links"]["events"] == "/v1/jobs/events"
            record = client.app.state.job_manager.wait(ref["job_id"])
        assert record.params == {
            "server_url": None,
            "genome_name": "rCRSd",
            "genome_digest": None,
            "asset_group_name": "fasta",
            "asset_name": None,
            "force": False,
        }

    def test_duplicate_submission_coalesces_with_202(self):
        rgc = web_stub_rgc()
        release = threading.Event()
        asset = rgc.pull.return_value
        rgc.pull.side_effect = lambda **kwargs: (release.wait(10), asset)[1]
        body = {"asset_group": "fasta", "genome": "rCRSd"}
        with make_local_client(rgc) as client:
            first = act(client, "POST", "/v1/actions/pull", json=body).json()
            second_response = act(client, "POST", "/v1/actions/pull", json=body)
            release.set()
            client.app.state.job_manager.wait(first["job_id"])
        assert second_response.status_code == 202  # there is no 409 job_in_progress
        second = second_response.json()
        assert second["job_id"] == first["job_id"]
        assert second["duplicate"] is True

    @pytest.mark.parametrize("server_url", [None, "http://mock-server"])
    def test_facade_call_contract(self, server_url):
        """server_url -> force_server_urls, and the invisible half: the runner
        must pass force_large=True and an explicit confirmer, or a pull can
        block a worker thread on a prompt nobody sees."""
        rgc = web_stub_rgc()
        body = {"asset_group": "fasta", "genome": "rCRSd"}
        if server_url:
            body["server_url"] = server_url
        with make_local_client(rgc) as client:
            ref = act(client, "POST", "/v1/actions/pull", json=body).json()
            client.app.state.job_manager.wait(ref["job_id"])
        kwargs = rgc.pull.call_args.kwargs
        assert kwargs["force_server_urls"] == ([server_url] if server_url else None)
        assert kwargs["force_large"] is True
        assert kwargs["confirm"] is not None

    @pytest.mark.parametrize(
        "body",
        [
            {"asset_group": "fasta"},  # neither genome nor genome_digest
            {"asset_group": "fasta", "genome": "g", "genome_digest": "d"},  # both
            {"asset_group": "fasta", "genome": "g", "bogus_field": 1},  # extra=forbid
        ],
    )
    def test_invalid_bodies_are_422_with_envelope(self, body):
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = act(client, "POST", "/v1/actions/pull", json=body)
        assert response.status_code == 422
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "validation_error"


class TestBuildEndpoint:
    def test_build_returns_202_and_converts_params(self):
        from refgenie.models import BuildParams

        rgc = web_stub_rgc()
        body = {
            "recipe": "bwa_index",
            "genome": "rCRSd",
            "asset_group": "bwa_index",
            "pull_parents": True,
            "params": {"files": {"data": "/tmp/data.txt"}, "params": {"cores": 4}},
        }
        with make_local_client(rgc) as client:
            response = act(client, "POST", "/v1/actions/build", json=body)
            assert response.status_code == 202
            ref = response.json()
            assert ref["kind"] == "build"
            record = client.app.state.job_manager.wait(ref["job_id"])
        assert record.status == "succeeded"
        kwargs = rgc.build_asset.call_args.kwargs
        assert kwargs["pull_parents"] is True
        assert isinstance(kwargs["params"], BuildParams)
        assert kwargs["params"].params == {"cores": 4}
        assert str(kwargs["params"].files["data"]) == "/tmp/data.txt"

    def test_preflight_reaches_the_facade_and_submits_no_job(self):
        """The router must use ``Refgenie.preflight_build``, never an
        ``AssetBuilder._``-private -- patching the facade proves the path."""
        rgc = web_stub_rgc()
        rgc.preflight_build.return_value = {
            "ok": True,
            "errors": [],
            "resolved": {"genome_digest": "genomedigest123", "asset_name": "default"},
        }
        with make_local_client(rgc) as client:
            response = act(
                client, "POST", "/v1/actions/build/preflight",
                json={"recipe": "fasta", "genome": "rCRSd", "asset_group": "fasta"},
            )
            jobs = client.get("/v1/jobs").json()["items"]
        assert response.status_code == 200
        payload = response.json()
        assert payload["ok"] is True
        assert payload["resolved"]["genome_digest"] == "genomedigest123"
        rgc.preflight_build.assert_called_once()
        assert rgc.preflight_build.call_args.kwargs["recipe_name"] == "fasta"
        assert jobs == []

    def test_preflight_reports_field_scoped_errors_with_200(self):
        rgc = web_stub_rgc()
        rgc.preflight_build.return_value = {
            "ok": False,
            "errors": [
                {"field": "params", "code": "missing_build_input", "message": "Missing 'data'."}
            ],
            "resolved": {},
        }
        with make_local_client(rgc) as client:
            response = act(
                client, "POST", "/v1/actions/build/preflight",
                json={"recipe": "needsfile", "genome": "rCRSd", "asset_group": "grp"},
            )
        assert response.status_code == 200  # a preflight that found problems succeeded
        payload = response.json()
        assert payload["ok"] is False
        assert payload["errors"][0]["field"] == "params"
        assert payload["errors"][0]["code"] == "missing_build_input"


class TestGenomeInitEndpoint:
    def test_genome_init_returns_202_and_runs_initialize_and_build(self):
        rgc = web_stub_rgc()
        body = {"fasta": "/tmp/genome.fa", "aliases": ["mygenome"], "species": "H. testens"}
        with make_local_client(rgc) as client:
            response = act(client, "POST", "/v1/actions/genomes", json=body)
            assert response.status_code == 202
            ref = response.json()
            assert ref["kind"] == "genome_init"
            record = client.app.state.job_manager.wait(ref["job_id"])
        assert record.status == "succeeded"
        assert record.result.genome_digest == "genomedigest123"
        kwargs = rgc.initialize_and_build.call_args.kwargs
        assert kwargs["genome_names"] == ["mygenome"]
        assert kwargs["species_name"] == "H. testens"


class TestSyncActions:
    """Real Refgenie, real filesystem, no mocks."""

    pytestmark = pytest.mark.component

    @pytest.fixture()
    def world(self, refgenie_built):
        with make_local_client(refgenie_built) as client:
            yield client, refgenie_built

    def _built_digest(self, rg):
        genome_digest = rg.alias.resolve("rCRSd")
        return genome_digest, rg.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="test"
        ).digest

    def test_delete_asset_removes_it(self, world):
        """Closes the untested-delete gap from the old dash pages."""
        client, rg = world
        genome_digest, asset_digest = self._built_digest(rg)
        response = act(client, "DELETE", f"/v1/actions/assets/{asset_digest}")
        assert response.status_code == 200
        payload = response.json()
        assert payload["ok"] is True
        assert payload["data"]["registry_path"]
        assert not rg.asset.exists(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="test"
        )

    def test_delete_asset_unknown_digest_is_404(self, world):
        client, _ = world
        response = act(client, "DELETE", "/v1/actions/assets/no_such_digest")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "asset_not_found"

    def test_delete_asset_with_children_is_409_naming_them(self, world, fixtures_path):
        client, rg = world
        genome_digest, parent_digest = self._built_digest(rg)
        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "demo.fa",
            alias_names=["demo"],
            description="demo genome",
        )
        rg.build_asset(
            recipe_name="fasta", genome_name="demo", asset_group_name="fasta", asset_name="test"
        )
        child = rg.asset.get(
            genome_digest=rg.alias.resolve("demo"), asset_group_name="fasta", asset_name="test"
        )
        rg.asset._asset_relations.set_children(genome_digest, "fasta", "test", [child.digest])

        response = act(client, "DELETE", f"/v1/actions/assets/{parent_digest}")
        assert response.status_code == 409
        error = response.json()["error"]
        assert error["code"] == "conflict"
        assert f"{rg.alias.resolve('demo')}/fasta" in error["message"]

    def test_delete_genome_by_alias_and_by_digest(self, world, fixtures_path):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        response = act(client, "DELETE", "/v1/actions/genomes/rCRSd")
        assert response.status_code == 200
        assert not rg.genome.exists(digest)

        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "demo.fa",
            alias_names=["demo"],
            description="demo genome",
        )
        demo_digest = rg.alias.resolve("demo")
        response = act(client, "DELETE", f"/v1/actions/genomes/{demo_digest}")
        assert response.status_code == 200
        assert not rg.genome.exists(demo_digest)

    def test_delete_genome_unknown_is_404(self, world):
        client, _ = world
        response = act(client, "DELETE", "/v1/actions/genomes/never_heard_of_it")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "genome_not_found"

    def test_alias_set_resolves_afterwards(self, world):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        response = act(
            client, "POST", "/v1/actions/aliases",
            json={"alias": "rCRSd2", "genome_digest": digest},
        )
        assert response.status_code == 200
        assert rg.alias.resolve("rCRSd2") == digest

    def test_alias_set_for_unknown_digest_is_404_not_a_phantom_genome(self, world):
        """Regression guard for the phantom-genome bug: without the existence
        guard, set_genome_alias silently creates a genome row for any string."""
        client, rg = world
        ghost = "0000000000_no_such_genome"
        response = act(
            client, "POST", "/v1/actions/aliases",
            json={"alias": "ghost", "genome_digest": ghost},
        )
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "genome_not_found"
        assert not rg.genome.exists(ghost)

    def test_alias_remove(self, world):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        rg.set_genome_alias(alias_name="doomed", genome_digest=digest)
        response = act(client, "DELETE", "/v1/actions/aliases/doomed")
        assert response.status_code == 200
        response = act(client, "DELETE", "/v1/actions/aliases/doomed")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "alias_not_found"

    def test_subscribe_reset_and_unsubscribe(self, world):
        client, rg = world
        response = act(
            client, "POST", "/v1/actions/subscriptions",
            json={"server_urls": ["http://a.example"]},
        )
        assert response.status_code == 200
        assert "http://a.example" in response.json()["data"]["subscriptions"]

        response = act(
            client, "POST", "/v1/actions/subscriptions",
            json={"server_urls": ["http://b.example"], "reset": True},
        )
        assert response.json()["data"]["subscriptions"] == ["http://b.example"]

        response = act(
            client, "DELETE", "/v1/actions/subscriptions",
            json={"server_urls": ["http://b.example"]},
        )
        assert response.status_code == 200
        assert "http://b.example" not in response.json()["data"]["subscriptions"]
        assert list(rg.configuration.get_server_subscriptions()) == []

    def test_set_default_asset(self, world):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        response = act(
            client, "POST", "/v1/actions/assets/default",
            json={"genome_digest": digest, "asset_group": "fasta", "asset": "test"},
        )
        assert response.status_code == 200
        assert rg.asset.get_default("fasta", genome_digest=digest) == "test"

    def test_set_default_asset_unknown_name_is_404(self, world):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        response = act(
            client, "POST", "/v1/actions/assets/default",
            json={"genome_digest": digest, "asset_group": "fasta", "asset": "nope"},
        )
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "asset_not_found"


class TestPreflightReal:
    """``Refgenie.preflight_build`` over a real world (no HTTP mocking)."""

    pytestmark = pytest.mark.component

    def test_valid_build_preflights_ok(self, refgenie_fs):
        with make_local_client(refgenie_fs) as client:
            response = act(
                client, "POST", "/v1/actions/build/preflight",
                json={"recipe": "fasta", "genome": "rCRSd", "asset_group": "fasta"},
            )
        assert response.status_code == 200
        payload = response.json()
        assert payload["ok"] is True, payload
        assert payload["errors"] == []
        assert payload["resolved"]["genome_digest"] == refgenie_fs.alias.resolve("rCRSd")
        assert payload["resolved"]["asset_name"] == "default"

    def test_unknown_genome_and_recipe_are_field_scoped(self, refgenie_fs):
        with make_local_client(refgenie_fs) as client:
            response = act(
                client, "POST", "/v1/actions/build/preflight",
                json={"recipe": "no_such_recipe", "genome": "no_such_genome", "asset_group": "x"},
            )
        payload = response.json()
        assert response.status_code == 200
        assert payload["ok"] is False
        by_field = {error["field"]: error["code"] for error in payload["errors"]}
        assert by_field["genome"] == "genome_not_found"
        assert by_field["recipe"] == "recipe_not_found"

    def test_missing_required_file_is_field_scoped(self, refgenie_fs, tmp_path):
        recipe_yaml = tmp_path / "needsfile_recipe.yaml"
        recipe_yaml.write_text(
            "\n".join(
                [
                    "name: needsfile",
                    "version: 0.1.0",
                    "output_asset_class: fasta",
                    "description: test recipe requiring an input file",
                    "input_files:",
                    "  data:",
                    "    description: required data file",
                    "input_params: null",
                    "input_assets: null",
                    "docker_image: null",
                    "command_templates:",
                    "  - cp {{values.files.data}} {{values.output_folder}}/",
                    'default_asset: "default"',
                ]
            )
        )
        refgenie_fs.recipe.add(recipe_yaml)
        with make_local_client(refgenie_fs) as client:
            response = act(
                client, "POST", "/v1/actions/build/preflight",
                json={"recipe": "needsfile", "genome": "rCRSd", "asset_group": "grp"},
            )
        payload = response.json()
        assert response.status_code == 200
        assert payload["ok"] is False
        assert any(
            e["field"] == "params" and e["code"] == "missing_build_input" and "data" in e["message"]
            for e in payload["errors"]
        ), payload["errors"]


class TestActionContract:
    """The surface is exactly ``EXPECTED`` -- a silently added or dropped
    endpoint is a build failure, the lesson of the rotted v1 tree."""

    EXPECTED = {
        ("POST", "/v1/actions/pull"),
        ("POST", "/v1/actions/build"),
        ("POST", "/v1/actions/build/preflight"),
        ("POST", "/v1/actions/genomes"),
        ("DELETE", "/v1/actions/assets/{asset_digest}"),
        ("DELETE", "/v1/actions/genomes/{genome_ref}"),
        ("POST", "/v1/actions/aliases"),
        ("DELETE", "/v1/actions/aliases/{alias_name}"),
        ("POST", "/v1/actions/subscriptions"),
        ("DELETE", "/v1/actions/subscriptions"),
        ("POST", "/v1/actions/assets/default"),
    }

    def _actions_routes(self):
        app = make_local_app(stub_rgc())
        return {
            (method, path)
            for method, path in _route_keys(app)
            if path.startswith("/v1/actions") and method != "HEAD"
        }

    def test_route_inventory_matches_the_contract(self):
        assert self._actions_routes() == self.EXPECTED


# ===========================================================================
# The localhost bridge: /ping, bridge modes, the LNA preflight
# ===========================================================================

#: Shared capability vocabulary -- the same fifteen keys /service-info emits.
CAPABILITY_KEY_SET = {
    "pull", "build", "delete", "aliases_write", "subscriptions", "recipes_write",
    "asset_classes_write", "remote_browse", "genome_init", "jobs", "jobs_cancel",
    "downloads", "archives", "seqcol", "drs",
}

REQUIRED_PING_FIELDS = {
    "service", "bridge_version", "mode", "refgenie_version", "api_version",
    "instance_id", "instance_label", "bridge_mode", "action_header",
    "capabilities", "bridge",
}


class TestPingContract:
    @pytest.mark.parametrize(
        # bridge_mode: "read" is the local default; the public API has no bridge.
        "mode,bridge_mode,pull,archives",
        [("local", "read", True, False), ("server", "off", False, True)],
    )
    def test_ping_shape(self, mode, bridge_mode, pull, archives):
        make_client = make_local_client if mode == "local" else make_server_client
        rgc = web_stub_rgc() if mode == "local" else stub_rgc()
        with make_client(rgc) as client:
            response = client.get("/ping")
        assert response.status_code == 200
        assert response.headers["cache-control"] == "no-store"
        payload = response.json()
        assert REQUIRED_PING_FIELDS <= set(payload)
        assert payload["service"] == "refgenie"
        assert isinstance(payload["bridge_version"], int)
        assert payload["action_header"] == "X-Refgenie-Action"
        assert payload["mode"] == mode
        assert payload["bridge_mode"] == bridge_mode
        assert payload["bridge"] == {"actions_cross_origin": False}
        assert set(payload["capabilities"]) == CAPABILITY_KEY_SET
        assert payload["capabilities"]["pull"] is pull
        assert payload["capabilities"]["archives"] is archives

    def test_ping_capabilities_match_service_info(self):
        """One shared vocabulary: /ping must emit exactly what /service-info
        emits, with no bridge-specific renames."""
        with make_local_client(web_stub_rgc()) as client:
            ping = client.get("/ping").json()
            service_info = client.get("/service-info").json()
        assert ping["capabilities"] == service_info["refgenie"]["capabilities"]

    def test_ping_reports_full_mode(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "full")
        with make_local_client(web_stub_rgc()) as client:
            payload = client.get("/ping").json()
        assert payload["bridge_mode"] == "full"
        assert payload["bridge"] == {"actions_cross_origin": True}

    def test_ping_still_answers_under_bridge_off(self, monkeypatch):
        """Same-origin callers (the local SPA) use /ping too; off only means no
        cross-origin caller can *read* it."""
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        with make_local_client(web_stub_rgc()) as client:
            payload = client.get("/ping").json()
        assert payload["bridge_mode"] == "off"

    def test_ping_omits_paths_by_default(self, tmp_path):
        rgc = web_stub_rgc()
        rgc.genome_folder = tmp_path / "genomes"
        with make_local_client(rgc) as client:
            payload = client.get("/ping").json()
        assert payload["instance_label"] == "local refgenie"
        assert str(tmp_path) not in str(payload)

    def test_ping_includes_paths_when_opted_in(self, monkeypatch, tmp_path):
        monkeypatch.setenv("REFGENIE_BRIDGE_EXPOSE_PATHS", "true")
        rgc = web_stub_rgc()
        rgc.genome_folder = tmp_path / "genomes"
        with make_local_client(rgc) as client:
            payload = client.get("/ping").json()
        assert payload["instance_label"] == str(tmp_path / "genomes")

    def test_instance_id_is_stable_across_app_constructions(self, monkeypatch, tmp_path):
        monkeypatch.setenv("REFGENIE_HOME_PATH", str(tmp_path))
        ids = []
        for _ in range(2):
            with make_local_client(web_stub_rgc()) as client:
                ids.append(client.get("/ping").json()["instance_id"])
        assert ids[0] == ids[1]
        assert (tmp_path / "instance_id").read_text().strip() == ids[0]


class TestBridgeCors:
    def test_preflight_from_bridge_origin_is_granted_by_default(self):
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(client, "/ping", PUBLIC_ORIGIN)
        assert response.status_code == 200
        assert response.headers["access-control-allow-origin"] == PUBLIC_ORIGIN

    def test_preflight_from_unlisted_origin_gets_no_grant(self):
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(client, "/ping", EVIL_ORIGIN)
        assert "access-control-allow-origin" not in response.headers

    def test_actual_read_carries_cors_headers_for_bridge_origin(self):
        """Job-status polling is a read; under ``read`` mode the bridge origin
        must be able to read /ping, /v4 and /v1/jobs responses."""
        with make_local_client(web_stub_rgc()) as client:
            for path in ("/ping", "/v1/jobs"):
                response = client.get(path, headers={"Origin": PUBLIC_ORIGIN})
                assert response.status_code == 200, path
                assert response.headers.get("access-control-allow-origin") == PUBLIC_ORIGIN, path

    def test_bridge_off_produces_no_cors_for_the_public_origin(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        with make_local_client(web_stub_rgc()) as client:
            options = preflight(client, "/ping", PUBLIC_ORIGIN)
            actual = client.get("/ping", headers={"Origin": PUBLIC_ORIGIN})
            same_origin = client.get("/ping")
        assert "access-control-allow-origin" not in options.headers
        assert "access-control-allow-origin" not in actual.headers
        assert same_origin.status_code == 200

    def test_origin_regex_escape_hatch(self, monkeypatch):
        monkeypatch.setenv(
            "REFGENIE_BRIDGE_ORIGIN_REGEX", r"https://[a-z0-9-]+\.refgenie-ui\.pages\.dev"
        )
        with make_local_client(web_stub_rgc()) as client:
            granted = client.get(
                "/ping", headers={"Origin": "https://deadbeef.refgenie-ui.pages.dev"}
            )
            refused = client.get("/ping", headers={"Origin": EVIL_ORIGIN})
        assert (
            granted.headers.get("access-control-allow-origin")
            == "https://deadbeef.refgenie-ui.pages.dev"
        )
        assert "access-control-allow-origin" not in refused.headers


class TestLocalNetworkAccessPreflight:
    @pytest.mark.parametrize(
        "request_header,response_header",
        [
            ("Access-Control-Request-Private-Network", "access-control-allow-private-network"),
            ("Access-Control-Request-Local-Network", "access-control-allow-local-network"),
        ],
    )
    def test_lna_grant_for_allowlisted_origin(self, request_header, response_header):
        """Both header spellings (PNA-era and LNA-era) are feature-detected. The
        status assertion is load-bearing: Starlette 400s any preflight carrying
        Access-Control-Request-Private-Network unless CORSMiddleware is built
        with ``allow_private_network=True``."""
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(client, "/ping", PUBLIC_ORIGIN, **{request_header: "true"})
        assert response.status_code == 200
        assert response.headers.get("access-control-allow-origin") == PUBLIC_ORIGIN
        assert response.headers.get(response_header) == "true"

    def test_lna_grant_withheld_for_unlisted_origin(self):
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(
                client, "/ping", EVIL_ORIGIN,
                **{"Access-Control-Request-Private-Network": "true"},
            )
        assert response.status_code == 400
        assert "access-control-allow-origin" not in response.headers

    def test_lna_middleware_absent_under_bridge_off(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(
                client, "/ping", PUBLIC_ORIGIN,
                **{"Access-Control-Request-Private-Network": "true"},
            )
        assert "access-control-allow-private-network" not in response.headers


PULL_BODY = {"asset_group": "fasta", "genome": "rCRSd"}


class TestCrossOriginActionPolicy:
    @pytest.mark.parametrize("bridge_mode", ["off", "read", "full"])
    def test_pull_without_action_header_is_403_in_every_mode(self, monkeypatch, bridge_mode):
        """The bridge must never weaken the anti-CSRF header requirement."""
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", bridge_mode)
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = client.post("/v1/actions/pull", json=PULL_BODY)
        assert response.status_code == 403
        assert response.json()["error"]["code"] == "missing_action_header"

    def test_cross_origin_pull_is_403_under_read_with_remedy(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "read")
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = act(
                client, "POST", "/v1/actions/pull", json=PULL_BODY,
                headers={"Origin": PUBLIC_ORIGIN},
            )
        assert response.status_code == 403
        payload = response.json()
        assert payload["error"]["code"] == "forbidden_origin"
        assert "refgenie dash --bridge full" in payload["error"]["message"]

    def test_cross_origin_pull_is_accepted_under_full(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "full")
        with make_local_client(web_stub_rgc()) as client:
            response = act(
                client, "POST", "/v1/actions/pull", json=PULL_BODY,
                headers={"Origin": PUBLIC_ORIGIN},
            )
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])

    def test_same_origin_pull_is_unaffected_by_bridge_mode(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        with make_local_client(web_stub_rgc()) as client:
            response = act(
                client, "POST", "/v1/actions/pull", json=PULL_BODY,
                headers={"Origin": "http://localhost"},  # == the request's own origin
            )
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])

    @pytest.mark.parametrize("bridge_mode", ["off", "read", "full"])
    def test_cross_origin_delete_is_403_in_every_mode(self, monkeypatch, bridge_mode):
        """Destructive verbs are never on the cross-origin surface."""
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", bridge_mode)
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = act(
                client, "DELETE", "/v1/actions/assets/deadbeef", headers={"Origin": PUBLIC_ORIGIN}
            )
        assert response.status_code == 403
        assert response.json()["error"]["code"] == "forbidden_origin"

    @pytest.mark.parametrize("bridge_mode", ["read", "full"])
    def test_unlisted_origin_actions_are_403(self, monkeypatch, bridge_mode):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", bridge_mode)
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = act(
                client, "POST", "/v1/actions/pull", json=PULL_BODY,
                headers={"Origin": EVIL_ORIGIN},
            )
        assert response.status_code == 403
        assert response.json()["error"]["code"] == "forbidden_origin"

    def test_dev_origin_keeps_full_action_access(self):
        """The Vite dev server serves the *local* SPA; it is trusted like
        same-origin and is not subject to the bridge's pull-only policy."""
        with make_local_client(web_stub_rgc()) as client:
            response = act(
                client, "POST", "/v1/actions/pull", json=PULL_BODY,
                headers={"Origin": "http://localhost:5173"},
            )
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])


class TestModeSurfaces:
    def test_foreign_host_is_accepted_in_server_mode(self):
        """The public API is proxied under real hostnames; the DNS-rebinding
        guard must not exist there."""
        with make_server_client(stub_rgc()) as client:
            response = client.get("/ping", headers={"host": "evil.example.com"})
        assert response.status_code == 200


# ===========================================================================
# One factory, two modes: create_app(mode="server"|"local")
# ===========================================================================


@pytest.fixture
def rgc(tmp_path, fixtures_path):
    """A catalog holding one genome, shared by the app-modes and SPA tests."""
    return make_server_rgc(tmp_path, fixtures_path, genomes=[("web_digest", ["rCRSd"])])


@pytest.fixture
def server_client(rgc):
    with TestClient(make_server_app(rgc), raise_server_exceptions=True) as client:
        yield client


@pytest.fixture
def local_client(rgc):
    # base_url: the local app's Host-header guard admits loopback names only,
    # so the default "testserver" host would 421 everything.
    with TestClient(
        make_local_app(rgc), base_url="http://localhost", raise_server_exceptions=True
    ) as client:
        yield client


class TestAppConstruction:
    """Both modes build, and neither reaches for the caller's real config."""

    @pytest.mark.parametrize("mode", [APP_MODE_SERVER, APP_MODE_LOCAL])
    def test_construction_does_not_touch_the_real_home(self, mode, rgc, tmp_path, monkeypatch):
        monkeypatch.setenv("HOME", str(tmp_path / "fake_home"))
        app = create_app(mode=mode, refgenie_instance=rgc)
        assert app.state.mode == mode
        assert not (tmp_path / "fake_home").exists()

    def test_unknown_mode_is_rejected(self, rgc):
        with pytest.raises(ValueError):
            create_app(mode="dashboard", refgenie_instance=rgc)


class TestSharedSurface:
    """The JSON API is the same in both modes, at the same prefix."""

    @pytest.mark.parametrize("client_name", ["server_client", "local_client"])
    def test_v4_genomes_lists_the_seeded_genome(self, client_name, request):
        client = request.getfixturevalue(client_name)
        response = client.get("/v4/genomes")
        assert response.status_code == 200
        assert [g["digest"] for g in response.json()["items"]] == ["web_digest"]


class TestServiceInfoBootstrap:
    """/service-info is how one SPA bundle serves both modes."""

    def test_reports_its_mode(self, server_client, local_client):
        assert server_client.get("/service-info").json()["refgenie"]["mode"] == "server"
        assert local_client.get("/service-info").json()["refgenie"]["mode"] == "local"

    def test_ga4gh_fields_stay_at_the_top_level(self, server_client):
        data = server_client.get("/service-info").json()
        assert data["id"] == "org.refgenie.api"
        assert data["type"]["group"] == "org.refgenie"
        assert "mode" not in data

    def test_local_mode_omits_the_seqcol_block(self, local_client, server_client):
        assert "seqcol" not in local_client.get("/service-info").json()
        assert "seqcol" in server_client.get("/service-info").json()

    def test_capability_key_set_is_identical_across_modes(self, server_client, local_client):
        """One vocabulary, one shape: the SPA reads the same keys in both modes."""
        server_caps = server_client.get("/service-info").json()["refgenie"]["capabilities"]
        local_caps = local_client.get("/service-info").json()["refgenie"]["capabilities"]
        assert set(server_caps) == set(local_caps) == CAPABILITY_KEY_SET
        assert all(isinstance(v, bool) for v in server_caps.values())

    #: Commands are local-only, bulk data is server-only, and the two
    #: definition writes are off everywhere -- the modes are mirror images.
    COMMAND_KEYS = ("pull", "build", "delete", "jobs", "remote_browse")
    DATA_KEYS = ("downloads", "archives", "seqcol", "drs")

    @pytest.mark.parametrize(
        "mode,commands_on", [(APP_MODE_SERVER, False), (APP_MODE_LOCAL, True)]
    )
    def test_capability_matrix(self, mode, commands_on):
        caps = _capabilities(mode)
        assert all(caps[key] is commands_on for key in self.COMMAND_KEYS)
        assert all(caps[key] is not commands_on for key in self.DATA_KEYS)
        assert not caps["recipes_write"] and not caps["asset_classes_write"]


#: (path, status in server mode, status in local mode).
_MODE_ISOLATION = [
    ("/v1/remote/genomes", 404, 200),
    ("/v1/remote/servers", 404, 200),
    ("/v1/jobs", 404, 200),
    ("/ga4gh/drs/service-info", 200, 404),
    ("/v4/ga4gh/drs/service-info", 200, 404),
    ("/seqcol/service-info", 200, 404),
    ("/data_channel/", 200, 404),
    ("/v4/archives", 200, 404),
]


@pytest.mark.parametrize("path,server_status,local_status", _MODE_ISOLATION)
def test_mode_isolation(server_client, local_client, path, server_status, local_status):
    # /v1/actions/* is absent from this table on purpose: the actions router has
    # no GET routes and the SPA catch-all answers GET for unmatched paths in
    # both modes, so a GET probe cannot tell the modes apart. Its mode isolation
    # is asserted by route inventory (TestActionContract,
    # TestActionSecurity.test_server_mode_has_no_actions_routes) and the
    # combined server-mode surface check below.
    assert server_client.get(path).status_code == server_status, f"server mode: {path}"
    assert local_client.get(path).status_code == local_status, f"local mode: {path}"


class TestRootNamespaceIsReservedForTheSpa:
    """The JSON API lives at /v4 and /v1 -- never at the root. The root prefix
    used to carry a second copy of the whole shared router in both apps."""

    @pytest.mark.parametrize("client_name", ["server_client", "local_client"])
    def test_root_genomes_is_not_the_api(self, client_name, request):
        """GET /genomes is an SPA page, not a JSON listing, in both modes."""
        client = request.getfixturevalue(client_name)
        response = client.get("/genomes")
        assert response.headers["content-type"].startswith("text/html")

    @pytest.mark.parametrize("client_name", ["server_client", "local_client"])
    def test_unmatched_api_path_is_a_json_404(self, client_name, request):
        """An API typo must not come back as a 200 HTML document -- and the body
        is the standard envelope, because the SPA catch-all builds this response
        itself, bypassing the handlers."""
        client = request.getfixturevalue(client_name)
        response = client.get("/v4/does_not_exist")
        assert response.status_code == 404
        assert response.headers["content-type"].startswith("application/json")
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "not_found"


class TestRouteHygiene:
    """Guards over the local app, matching the ones test_server.py runs."""

    def test_local_app_has_no_duplicate_routes(self, rgc):
        _assert_no_duplicate_routes(make_local_app(rgc), "local app")

    def test_local_app_has_unique_operation_ids(self, rgc):
        _assert_unique_operation_ids(make_local_app(rgc), "local app")

    @pytest.mark.parametrize("segment", SPA_CLIENT_ROUTES)
    def test_spa_client_routes_do_not_collide_with_the_api(self, rgc, segment):
        """No API route may shadow an SPA page, in either mode."""
        for app, label in ((make_server_app(rgc), "server"), (make_local_app(rgc), "local")):
            assert ("GET", f"/{segment}") not in set(_route_keys(app)), f"{label}: /{segment}"


#: The frontend route table. Present in a source checkout; absent from a wheel.
_ROUTES_TSX = Path(__file__).parent.parent / "frontend" / "src" / "app" / "routes.tsx"


@pytest.mark.skipif(not _ROUTES_TSX.is_file(), reason="frontend source not present")
def test_spa_route_mirror():
    """``SPA_CLIENT_ROUTES`` (Python) and the frontend route table must agree.
    Two dumb regex passes over ``routes.tsx``, on purpose -- no TS parser."""
    source = _ROUTES_TSX.read_text()

    const_block = re.search(r"export const SPA_CLIENT_ROUTES = \[(.*?)\]", source, re.DOTALL)
    assert const_block, "routes.tsx no longer exports SPA_CLIENT_ROUTES"
    ts_const = set(re.findall(r"'/([a-z0-9-]+)'", const_block.group(1)))
    assert ts_const == set(SPA_CLIENT_ROUTES), (
        "frontend/src/app/routes.tsx SPA_CLIENT_ROUTES drifted from "
        "refgenie/server/const.py:\n"
        f"  frontend only: {sorted(ts_const - set(SPA_CLIENT_ROUTES))}\n"
        f"  backend only:  {sorted(set(SPA_CLIENT_ROUTES) - ts_const)}"
    )

    registered = re.findall(r"path:\s*'([^']+)'", source)
    segments = {path.split("/")[0] for path in registered} - {"", "*"}
    assert segments == set(SPA_CLIENT_ROUTES), (
        "registered frontend routes drifted from SPA_CLIENT_ROUTES:\n"
        f"  registered only: {sorted(segments - set(SPA_CLIENT_ROUTES))}\n"
        f"  tuple only:      {sorted(set(SPA_CLIENT_ROUTES) - segments)}"
    )


# ===========================================================================
# SPA serving (refgenie.server.spa)
# ===========================================================================

#: Matches a hashed JS asset reference the way vite emits it.
_HASHED_JS_RE = re.compile(r"/_app/[A-Za-z0-9._-]+\.js")


@requires_web_assets
class TestSpaServing:
    """Everything here needs the frontend built into ``refgenie/server/webui/``."""

    @pytest.mark.parametrize("client_name", ["server_client", "local_client"])
    def test_spa_index_served(self, client_name, request):
        """Both modes serve the SPA shell at ``/`` -- there is no mode-specific root."""
        client = request.getfixturevalue(client_name)
        response = client.get("/")
        assert response.status_code == 200
        assert response.headers["content-type"].startswith("text/html")
        assert 'id="root"' in response.text

    def test_hashed_asset_served(self, server_client):
        """The hash is derived from the index, never hard-coded -- it rots on rebuild."""
        index = server_client.get("/").text
        match = _HASHED_JS_RE.search(index)
        assert match, f"no hashed JS asset referenced in index.html: {index!r}"
        response = server_client.get(match.group(0))
        assert response.status_code == 200
        assert response.content
        assert "javascript" in response.headers["content-type"]

    def test_asset_cache_headers(self, server_client):
        index_response = server_client.get("/")
        assert "no-cache" in index_response.headers["cache-control"]

        match = _HASHED_JS_RE.search(index_response.text)
        assert match
        asset_response = server_client.get(match.group(0))
        assert "immutable" in asset_response.headers["cache-control"]

    def test_client_route_falls_back_to_index(self, local_client):
        """Deep links are the first thing users hit and the first thing to break."""
        index = local_client.get("/").text
        deep_link = local_client.get("/genomes/abc123")
        assert deep_link.status_code == 200
        assert deep_link.text == index

    def test_unknown_asset_404s(self, server_client):
        """The registration-order regression test: a missing hashed asset must
        never fall back to index.html with a 200."""
        response = server_client.get("/_app/nope-deadbeef.js")
        assert response.status_code == 404
        assert "text/html" not in response.headers["content-type"]

    def test_service_info_bootstrap(self, server_client):
        """The SPA reads its bootstrap from /service-info (no injected global)."""
        data = server_client.get("/service-info").json()
        refgenie = data["refgenie"]
        assert refgenie["mode"] == "server"
        assert refgenie["api_base"] == "/v4"
        assert refgenie["root_path"] == ""
        assert "web_ui" in refgenie
        assert set(refgenie["capabilities"]) == set(CAPABILITY_KEYS)
        assert all(isinstance(v, bool) for v in refgenie["capabilities"].values())

    def test_root_path_rewrites_base_href(self, rgc):
        """One build, every deployment: root_path is a construction-time rewrite,
        not a rebuild. TestClient talks to the app post-reverse-proxy."""
        app = create_app(mode=APP_MODE_SERVER, root_path="/refgenie", refgenie_instance=rgc)
        with TestClient(app, raise_server_exceptions=True) as client:
            response = client.get("/")
            info = client.get("/service-info").json()
        assert '<base href="/refgenie/" />' in response.text
        assert info["refgenie"]["root_path"] == "/refgenie"

    def test_route_guard_single_spa_mount(self, rgc):
        """Exactly one catch-all, registered last -- Starlette matches routes in
        registration order, so anything after it would be dead."""
        app = create_app(mode=APP_MODE_SERVER, refgenie_instance=rgc)
        catch_alls = [
            r for r in app.router.routes if getattr(r, "path", None) == "/{full_path:path}"
        ]
        assert len(catch_alls) == 1
        assert app.router.routes[-1] is catch_alls[0]


def test_missing_assets_returns_503(rgc, monkeypatch):
    """Absence must be tolerated, never fatal. Not marked ``requires_web_assets``:
    this is the Node-free contributor path, exercised regardless of whether this
    checkout has a bundle. Patch the resolver itself so the "no bundle anywhere"
    case is exercised even on a machine that has built the frontend."""
    import refgenie.server.main as server_main

    monkeypatch.setattr(server_main, "resolve_web_dist", lambda explicit=None: None)

    app = create_app(mode=APP_MODE_LOCAL, refgenie_instance=rgc)
    with TestClient(app, base_url="http://localhost", raise_server_exceptions=True) as client:
        index_response = client.get("/")
        api_response = client.get("/v4/genomes")

    assert index_response.status_code == 503
    assert "npm --prefix frontend run build" in index_response.text
    assert api_response.status_code == 200


# ===========================================================================
# OpenAPI <-> TypeScript wire-type drift guard (frontend/src/types/api.ts)
# ===========================================================================

#: Property key sets, mirroring frontend/src/types/api.ts one interface at a time.
EXPECTED_PROPERTIES = {
    "GenomeResponse": {
        "digest", "aliases", "description", "asset_count", "species_name",
        "common_name", "taxon_id", "assembly_source", "assembly_accession",
    },
    "GenomeDetailResponse": {
        "digest", "description", "species_name", "common_name", "taxon_id",
        "assembly_source", "assembly_accession", "assembly_level", "remote_url",
        "taxon_uri", "fhr",
    },
    "AssetGroupPublic": {"id", "name", "description", "genome_digest", "asset_class_id"},
    "AssetResponse": {
        "digest", "name", "description", "recipe_id", "asset_group_id", "size",
        "serving_modes_override", "colocate", "serving_modes", "asset_class_name",
        "asset_group_name", "genome_digest", "names", "seek_keys", "is_default",
    },
    "SeekKeyResponse": {"name", "value", "description", "type", "size"},
    "AssetClassPublic": {"id", "name", "version", "description", "serving_modes"},
    "RecipePublic": {
        "id", "name", "version", "description", "output_asset_class_id",
        "command_templates", "input_params", "input_files", "input_assets",
        "docker_image", "custom_seek_keys", "default_asset", "inherent",
    },
    "StagedAssetPublic": {
        "asset_digest", "mode", "directory_contents", "build_commands",
        "download_count", "tarball_digest", "tarball_size",
    },
    "AliasPublic": {"name", "genome_digest"},
    "AliasResponse": {"alias", "digest", "source", "collection", "fhr"},
    "PaginationMeta": {"offset", "limit", "total"},
    "ArchiveRecord": {
        "digest", "asset_digest", "tarball_digest", "size", "directory_contents",
        "build_commands", "download_count",
    },
    "DatabaseSummaryResponse": {"genomes", "asset_groups", "assets"},
    "SpeciesStatistics": {"genomes", "asset_classes", "assets"},
}

#: Endpoints the SPA's service layer calls on the shared read API.
EXPECTED_SHARED_PATHS = {
    "/genomes", "/genomes/{digest}", "/asset_groups", "/asset_groups/{id}",
    "/assets", "/assets/{digest}", "/assets/{asset_digest}/files", "/asset_classes",
    "/asset_classes/{id}", "/recipes", "/recipes/{id}", "/configurations",
    "/configurations/{id}", "/staged_assets", "/staged_assets/{id}",
    "/relationships/{asset_digest}", "/aliases", "/aliases/{name}",
}

#: Server-only reads the SPA gates on capabilities.archives / .downloads.
EXPECTED_VERSION4_PATHS = {
    "/archives", "/archives/{asset_digest}/download",
    "/assets/{asset_digest}/files/{file_path}", "/species/summary", "/summary",
}


@pytest.fixture(scope="module")
def openapi_schema():
    app = FastAPI()
    app.include_router(shared_router.router, prefix="/v4")
    app.include_router(version4_router.router, prefix="/v4")
    return app.openapi()


@pytest.mark.parametrize("model_name", sorted(EXPECTED_PROPERTIES))
def test_wire_model_properties_match_frontend_types(openapi_schema, model_name):
    schemas = openapi_schema["components"]["schemas"]
    assert model_name in schemas, (
        f"{model_name} is no longer in the OpenAPI schema. "
        "frontend/src/types/api.ts declares it; update both."
    )
    actual = set(schemas[model_name].get("properties", {}))
    expected = EXPECTED_PROPERTIES[model_name]
    assert actual == expected, (
        f"{model_name} properties drifted from frontend/src/types/api.ts.\n"
        f"  added on the server:   {sorted(actual - expected)}\n"
        f"  missing on the server: {sorted(expected - actual)}\n"
        "Update frontend/src/types/api.ts (and whatever renders the field), "
        "then update EXPECTED_PROPERTIES here."
    )


@pytest.mark.parametrize(
    "expected,service_module",
    [
        (EXPECTED_SHARED_PATHS, "resources/*.ts"),
        (EXPECTED_VERSION4_PATHS, "resources/serverInfo.ts"),
    ],
    ids=["shared", "server_only"],
)
def test_paths_the_spa_calls_are_stable(openapi_schema, expected, service_module):
    paths = {p.removeprefix("/v4") for p in openapi_schema["paths"]}
    missing = expected - paths
    assert not missing, (
        f"Endpoints the SPA calls disappeared: {sorted(missing)}. "
        f"Update frontend/src/services/{service_module}."
    )


@pytest.mark.parametrize(
    "enum_name,members",
    [
        ("SeekKeyType", {"file", "directory", "prefix", "string", "json"}),
        ("SearchOperator", {"eq", "contains", "starts_with", "ends_with"}),
    ],
)
def test_enums_match_the_frontend_unions(openapi_schema, enum_name, members):
    assert set(openapi_schema["components"]["schemas"][enum_name]["enum"]) == members


# ===========================================================================
# Building from a worker thread; failing fast on a bad stage request
# ===========================================================================


class TestBuildOnWorkerThread:
    """``AssetBuilder.build`` registered a SIGINT handler unconditionally;
    ``signal.signal()`` raises off the main thread, so every threaded build died
    on its first line of real work. Both tests build real assets: component."""

    pytestmark = pytest.mark.component

    def test_build_succeeds_on_a_worker_thread(self, refgenie_fs):
        """The anchor test for the main-thread guard on ``signal.signal``."""
        box = {}

        def build():
            try:
                box["asset"] = refgenie_fs.build_asset(
                    recipe_name="fasta",
                    genome_name=GENOME,
                    asset_group_name=GROUP,
                    asset_name=ASSET,
                )
            except BaseException as exc:  # noqa: BLE001 - re-raised below
                box["error"] = exc

        thread = threading.Thread(target=build)
        thread.start()
        thread.join(timeout=120)
        assert not thread.is_alive(), "the build thread never finished"

        if "error" in box:
            raise box["error"]
        assert box["asset"] is not None
        assert box["asset"].digest

    def test_the_main_thread_still_registers_a_handler(self, refgenie_fs, monkeypatch):
        """The guard must not turn the CLI's Ctrl-C handling off."""
        import refgenie.managers.asset.builder as builder_module

        registered = []
        real_signal = builder_module.signal.signal

        def spy(sig, handler):
            registered.append(sig)
            return real_signal(sig, handler)

        monkeypatch.setattr(builder_module.signal, "signal", spy)
        refgenie_fs.build_asset(
            recipe_name="fasta", genome_name=GENOME, asset_group_name=GROUP, asset_name=ASSET
        )
        import signal as signal_module

        assert signal_module.SIGINT in registered


class TestStageValidationHappensFirst:
    """``stage=True`` with no stage folder must fail before anything is built.
    These use a test double, so they stay unit."""

    pytestmark = pytest.mark.unit

    @pytest.fixture
    def unstaged(self, refgenie_with_fasta, monkeypatch):
        """A Refgenie whose config names no stage folder."""
        monkeypatch.setattr(
            type(refgenie_with_fasta), "genome_stage_folder", property(lambda self: None)
        )
        return refgenie_with_fasta

    def test_it_raises_without_building(self, unstaged, monkeypatch):
        called = []
        monkeypatch.setattr(unstaged.asset, "build", lambda **kwargs: called.append(kwargs))

        with pytest.raises(ValueError, match="genome_stage_folder is not set"):
            unstaged.build_asset(
                recipe_name="fasta",
                genome_name=GENOME,
                asset_group_name=GROUP,
                asset_name=ASSET,
                stage=True,
            )
        assert called == [], "the build ran before the stage folder was validated"

    def test_no_staging_requested_is_unaffected(self, unstaged, monkeypatch):
        monkeypatch.setattr(unstaged.alias, "resolve", lambda name: "d" * 64)
        monkeypatch.setattr(unstaged.asset, "build", lambda **kwargs: None)

        assert (
            unstaged.build_asset(
                recipe_name="fasta", genome_name=GENOME, asset_group_name=GROUP, asset_name=ASSET
            )
            is None
        )


# ===========================================================================
# End-to-end: one real pull/build all the way through the job manager
# ===========================================================================


class TestJobsIntegration:
    """These prove the seams meet: the progress sink installed by the manager
    reaches ``download_with_progress`` four call layers down, the ``refgenie``
    logger's narration is attributed to the right job, and the asset really
    lands on disk. Real builds/pulls: component."""

    pytestmark = pytest.mark.component

    @pytest.fixture
    def manager(self, server_client_world):
        """A JobManager over the client half of a served server/client pair."""
        client_rg, _, _ = server_client_world
        manager = JobManager(client_rg)
        yield manager
        manager.shutdown(wait=False)

    def test_a_real_pull_reports_progress_logs_and_lands_the_asset(
        self, manager, server_client_world
    ):
        client_rg, _, url = server_client_world

        ref = manager.submit_pull(
            PullJobParams(server_url=url, genome_name=GENOME, asset_group_name=GROUP)
        )
        record = manager.wait(ref.job_id, timeout=120)

        assert record.status == JobStatus.SUCCEEDED, record.error
        assert record.result.registry_path
        assert record.result.asset_digest

        events, _, truncated = manager.events_since(0)
        assert truncated is False
        assert all(event.job_id == ref.job_id for event in events)

        byte_events = [e for e in events if e.type == "progress" and e.bytes_done is not None]
        assert byte_events, "the progress sink never reached download_with_progress"
        counts = [e.bytes_done for e in byte_events]
        assert counts == sorted(counts), "byte counts must not go backwards"

        phases = [e.phase for e in events if e.type == "progress" and e.phase]
        assert "download" in phases
        assert "verify" in phases
        assert "register" in phases
        assert set(phases) <= PHASE_VOCABULARY["pull"], "a phase outside the UI's vocabulary"

        log_lines = [e.line for e in events if e.type == "log"]
        assert log_lines, "no puller narration was captured"
        assert all(e.source == "refgenie" for e in events if e.type == "log")

        asset = client_rg.asset.get_by_digest(record.result.asset_digest)
        assert asset is not None
        assert (Path(client_rg.genome_folder) / asset.path).is_dir()

    def test_repulling_an_existing_asset_reports_asset_exists(self, manager, server_client_world):
        """The differentiated pull error survives the move to jobs: a real
        ``AssetExistsError`` must come back as ``asset_exists`` -- the one the
        UI answers with a "force" button -- not a generic failure."""
        _, _, url = server_client_world
        params = PullJobParams(server_url=url, genome_name=GENOME, asset_group_name=GROUP)

        first = manager.submit_pull(params)
        assert manager.wait(first.job_id, timeout=120).status == JobStatus.SUCCEEDED

        second = manager.submit_pull(params)
        assert second.job_id != first.job_id, "a finished job must not coalesce a new one"
        assert second.duplicate is False

        record = manager.wait(second.job_id, timeout=120)
        assert record.status == JobStatus.FAILED
        assert record.error.code == ErrorCode.ASSET_EXISTS


class TestLocalModeSmoke:
    """The assembled ``refgenie dash`` app, boot to job, over HTTP. The deep
    coverage lives above; this asserts that ``create_app(mode="local")`` wires
    the actions router, the job manager and the SSE surface together."""

    pytestmark = pytest.mark.component

    def test_job_lifecycle_through_the_real_runner(self, refgenie_fs):
        """Submit -> queued 202 -> real pull runner (no subscriptions, so it
        fails) -> terminal record with a classified code, readable over HTTP.
        A *successful* end-to-end pull is TestJobsIntegration's contract."""
        with make_local_client(refgenie_fs) as client:
            response = client.post(
                "/v1/actions/pull",
                json={"asset_group": "fasta", "genome": "rCRSd"},
                headers=ACTION_HEADERS,
            )
            assert response.status_code == 202
            ref = response.json()
            assert ref["status"] == "queued"

            record = client.app.state.job_manager.wait(ref["job_id"])
            assert record.status == "failed"
            assert record.error.code == "no_subscriptions"

            over_http = client.get(f"/v1/jobs/{ref['job_id']}").json()
        assert over_http["status"] == "failed"
        assert over_http["error"]["code"] == "no_subscriptions"


# ===========================================================================
# Regression guard: every CLI model field is consumed by its handler
# ===========================================================================

# Fields consumed indirectly, not via ``cmd.<field>`` in the handler body.
# Keep this short: growth here means a flag went unwired again.
_INDIRECT = {
    # PromptHandlingMixin reads these through ``self.<field>`` in
    # resolve_force_large() and its validators (refgenie/cli/commands/pull.py).
    "skip_large",
    "pull_large",
    "batch",
}


def _model_handler_pairs():
    """Yield (model, handler) for every leaf command model in the CLI."""
    pairs = []
    for model, handler in get_dispatch().items():
        group_dispatch = None
        for name, value in vars(inspect.getmodule(handler)).items():
            if name.endswith("_DISPATCH") and isinstance(value, dict):
                group_dispatch = value
        if model.__name__.endswith("Parser"):
            if group_dispatch:
                pairs.extend(group_dispatch.items())
        else:
            pairs.append((model, handler))
    # asset.handle_asset_group forwards its only leaf to listing.handle_list.
    pairs.append((cli_asset.AssetListNestedModel, cli_listing.handle_list))
    return pairs


_PAIRS = _model_handler_pairs()


@pytest.mark.parametrize(
    "model,handler",
    _PAIRS,
    ids=[f"{m.__name__}-{h.__name__}" for m, h in _PAIRS],
)
def test_every_model_field_is_consumed(model, handler):
    """A field with no handler reader renders a flag that is a silent no-op.
    The rule is wire-or-delete: if this fails, consume the field or remove it --
    do not extend ``_INDIRECT``."""
    source = inspect.getsource(handler)
    # A handler may delegate ``cmd`` to a shared helper; count what those
    # helpers read as consumed.
    for name, value in vars(cli_helpers).items():
        if inspect.isfunction(value) and re.search(rf"\b{name}\(", source):
            source += inspect.getsource(value)
    unwired = [
        field
        for field in model.model_fields
        if field not in _INDIRECT
        and not re.search(rf"cmd\.{field}\b", source)
        and not re.search(rf"""getattr\(\s*cmd\s*,\s*["']{field}["']""", source)
    ]
    assert not unwired, (
        f"{model.__name__} fields {unwired} are never read by "
        f"{handler.__module__}.{handler.__name__}. Wire the flag or delete the field."
    )
