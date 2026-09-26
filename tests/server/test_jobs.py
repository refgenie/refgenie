"""``JobManager`` and the jobs HTTP API (``refgenie/server/jobs/``).

Covers the manager's lifecycle, cancel, coalescing, history and ring
truncation against fake runners; the real pull/build runners' error mapping to
the shared error vocabulary, as job records; and the jobs HTTP API with its
one multiplexed SSE stream, replay, ``Last-Event-ID``, wire contract and log
capture.

Apps here are a bare ``FastAPI()`` with the real jobs router included -- on
purpose: the contract is that the jobs router carries no dependency on the app
factory. Its mode isolation *inside* ``create_app`` is asserted in
``tests/server/test_app.py``.
"""

import logging
import threading
import time
from unittest.mock import MagicMock

import pytest

from tests.helpers import (
    ACTION_HEADERS,
    BUILD_PARAMS,
    PULL_PARAMS,
    build_params,
    instant_runner,
    job_result_ok,
    jobs_app,
    jobs_client,
    make_gated_runner,
    make_slow_runner,
    open_gate,
    pull_params,
    read_sse,
    requires_dash,
    requires_server,
    submit,
    submit_and_wait,
    until_done,
    wait_for,
)

requires_server()
requires_dash()

from fastapi import FastAPI  # noqa: E402  (must follow the extras guard)
from fastapi.testclient import TestClient  # noqa: E402

from refgenie import progress  # noqa: E402
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
from refgenie.server.errors import ErrorCode, install_error_handlers  # noqa: E402
from refgenie.server.jobs import jobs_router  # noqa: E402
from refgenie.server.jobs.manager import CancelOutcome, JobManager  # noqa: E402
from refgenie.server.jobs.schemas import JobKind, JobStatus  # noqa: E402

#: The result every fake runner returns, as the typed model the manager stores.
OK = job_result_ok()


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
            (
                AssetExistsError("x"),
                JobStatus.FAILED,
                ErrorCode.ASSET_EXISTS,
                "Asset already exists. Use force to overwrite.",
            ),
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
            lambda: (
                manager.record(a.job_id).status == JobStatus.RUNNING
                and manager.record(b.job_id).status == JobStatus.RUNNING
            )
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
# Pull and build error mapping, as job records
# ===========================================================================

RUNNER_PULL = pull_params(server_url="http://mock-server", genome_name="test_genome", force=True)
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
            (
                AssetExistsError("test/fasta:default already exists"),
                JobStatus.FAILED,
                ErrorCode.ASSET_EXISTS,
            ),
            (
                NoArchiveError("No archive found for asset digest abc123"),
                JobStatus.FAILED,
                ErrorCode.NO_ARCHIVE,
            ),
            (PullFailedError("Server refused"), JobStatus.FAILED, ErrorCode.PULL_FAILED),
            (PullSkipped("Skipping pull of x"), JobStatus.CANCELLED, ErrorCode.PULL_SKIPPED),
            (RuntimeError("something unexpected"), JobStatus.FAILED, ErrorCode.INTERNAL_ERROR),
        ],
    )
    def test_exception_maps_to_status_and_code(self, run_job, exc, expected_status, expected_code):
        rgc = MagicMock()
        rgc.transfer.pull.side_effect = exc
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.status == expected_status
        assert record.error.code == expected_code
        assert record.error.detail, "a traceback belongs on the record"

    def test_asset_exists_message_is_actionable(self, run_job):
        rgc = MagicMock()
        rgc.transfer.pull.side_effect = AssetExistsError("test/fasta:default already exists")
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
        rgc.transfer.pull.side_effect = exc
        record = run_job("pull", RUNNER_PULL, rgc)
        assert fragment in record.error.message.lower()

    def test_returning_none_is_no_subscriptions(self, run_job):
        """``pull`` returns None when nothing is subscribed and the subscribe
        prompt was declined -- which, from a job, it always is."""
        rgc = MagicMock()
        rgc.transfer.pull.return_value = None
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.status == JobStatus.FAILED
        assert record.error.code == ErrorCode.NO_SUBSCRIPTIONS
        assert "subscribe" in record.error.message.lower()


class TestPullSuccess:
    def test_result_carries_what_the_ui_needs_to_refresh(self, run_job):
        rgc = MagicMock()
        rgc.transfer.pull.return_value = fake_asset()
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.status == JobStatus.SUCCEEDED
        assert record.result.asset_digest == "a" * 64
        assert record.result.registry_path == "dig/fasta:default"
        assert record.result.asset_group_name == "fasta"
        assert record.result.genome_digest == "g" * 64
        assert record.result.staged is False

    def test_the_root_is_called_non_interactively(self, run_job):
        """``force_large=True`` because a browser cannot answer a size prompt,
        and an explicit ``confirm`` because ``resolve_confirmer(None)`` falls
        through to a real ``rich`` prompt once the CLI has enabled prompts."""
        rgc = MagicMock()
        rgc.transfer.pull.return_value = fake_asset()
        run_job("pull", RUNNER_PULL, rgc)

        kwargs = rgc.transfer.pull.call_args.kwargs
        assert kwargs["force_large"] is True
        assert kwargs["confirm"] is not None
        assert kwargs["confirm"]("replace everything?") is False
        assert kwargs["force_server_urls"] == ["http://mock-server"]
        assert kwargs["genome"] == "test_genome"
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
        rgc.transfer.pull.return_value = asset
        record = run_job("pull", RUNNER_PULL, rgc)
        assert record.status == JobStatus.SUCCEEDED
        assert record.result.asset_digest == "a" * 64
        assert record.result.registry_path is None


class TestBuildErrors:
    def test_success(self, run_job, tmp_path):
        rgc = MagicMock()
        rgc.genome_folder = tmp_path
        rgc.build.run.return_value = fake_asset()
        record = run_job("build", RUNNER_BUILD, rgc)
        assert record.status == JobStatus.SUCCEEDED
        assert record.result.registry_path == "dig/fasta:default"
        assert record.result.asset_digest == "a" * 64

    @pytest.mark.parametrize("image, docker", [("docker.io/databio/refgenie", True), (None, False)])
    def test_build_runs_in_docker_when_the_recipe_declares_an_image(
        self, run_job, tmp_path, image, docker
    ):
        """Preflight names the asset inside the recipe's image; the build must
        run there too, or a container recipe preflights clean and then fails
        on a host without the tool. No image means a native build."""
        rgc = MagicMock()
        rgc.genome_folder = tmp_path
        rgc.recipe.get.return_value.docker_image = image
        rgc.build.run.return_value = fake_asset()
        record = run_job("build", RUNNER_BUILD, rgc)
        assert record.status == JobStatus.SUCCEEDED
        assert rgc.build.run.call_args.kwargs["docker"] is docker

    def test_none_return_is_build_failed(self, run_job, tmp_path):
        """``build.run`` returning None means the pipeline failed. The library
        answer is a printed message and exit 0; a job card needs a code."""
        rgc = MagicMock()
        rgc.genome_folder = tmp_path
        rgc.build.run.return_value = None
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

        rgc.build.run.side_effect = slow_build
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
    def test_build_exceptions_use_the_shared_vocabulary(
        self, run_job, tmp_path, exc, expected_code
    ):
        rgc = MagicMock()
        rgc.genome_folder = tmp_path
        rgc.build.run.side_effect = exc
        record = run_job("build", RUNNER_BUILD, rgc)
        assert record.status == JobStatus.FAILED
        assert record.error.code == expected_code


# ===========================================================================
# The jobs HTTP API, including the multiplexed SSE stream
#
# The app here is a bare FastAPI() with the real router included -- deliberate:
# the contract is that the jobs router carries no dependency on the app factory.
# Its mode isolation *inside* create_app is asserted in tests/server/test_app.py.
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
        assert empty == {"lines": [], "next_offset": 0, "has_more": False}

        log = tmp_path / "pipeline_log.md"
        log.write_text("one\ntwo\nthree\n")
        manager.get(job_id).log_file = str(log)

        first_page = client.get(f"/v1/jobs/{job_id}/log", params={"limit": 2}).json()
        assert first_page == {"lines": ["one", "two"], "next_offset": 2, "has_more": True}
        second_page = client.get(f"/v1/jobs/{job_id}/log", params={"offset": 2}).json()
        assert second_page == {"lines": ["three"], "next_offset": 3, "has_more": False}

    def test_log_negative_offset_is_422(self, client, manager):
        job_id = submit_and_wait(client, manager)
        response = client.get(f"/v1/jobs/{job_id}/log", params={"offset": -1})
        assert response.status_code == 422

    def test_log_paging_across_window_boundaries(self, client, manager, tmp_path):
        job_id = submit_and_wait(client, manager)
        lines = [f"l{i}" for i in range(10)]
        log = tmp_path / "pipeline_log.md"
        log.write_text("\n".join(lines) + "\n")
        manager.get(job_id).log_file = str(log)

        pages = []
        offsets = []
        has_mores = []
        offset = 0
        while True:
            body = client.get(
                f"/v1/jobs/{job_id}/log", params={"offset": offset, "limit": 3}
            ).json()
            pages.append(body["lines"])
            offsets.append(body["next_offset"])
            has_mores.append(body["has_more"])
            offset = body["next_offset"]
            if not body["has_more"]:
                break

        assert pages == [["l0", "l1", "l2"], ["l3", "l4", "l5"], ["l6", "l7", "l8"], ["l9"]]
        assert [line for page in pages for line in page] == lines
        assert offsets == [3, 6, 9, 10]
        assert has_mores == [True, True, True, False]

    def test_log_has_more_false_at_exact_end(self, client, manager, tmp_path):
        job_id = submit_and_wait(client, manager)
        lines = [f"l{i}" for i in range(6)]
        log = tmp_path / "pipeline_log.md"
        log.write_text("\n".join(lines) + "\n")
        manager.get(job_id).log_file = str(log)

        body = client.get(f"/v1/jobs/{job_id}/log", params={"offset": 3, "limit": 3}).json()
        assert body == {"lines": ["l3", "l4", "l5"], "next_offset": 6, "has_more": False}

        body = client.get(f"/v1/jobs/{job_id}/log", params={"offset": 6, "limit": 3}).json()
        assert body == {"lines": [], "next_offset": 6, "has_more": False}

        body = client.get(f"/v1/jobs/{job_id}/log", params={"offset": 100, "limit": 3}).json()
        assert body == {"lines": [], "next_offset": 100, "has_more": False}

    def test_log_file_without_trailing_newline(self, client, manager, tmp_path):
        job_id = submit_and_wait(client, manager)
        log = tmp_path / "pipeline_log.md"
        log.write_text("a\nb")
        manager.get(job_id).log_file = str(log)

        body = client.get(f"/v1/jobs/{job_id}/log").json()
        assert body["lines"] == ["a", "b"]

    def test_log_overlong_line_is_capped(self, client, manager, tmp_path):
        from refgenie.server.jobs.router import _LOG_LINE_MAX_CHARS, _LOG_LINE_TRUNCATED_SUFFIX

        job_id = submit_and_wait(client, manager)
        log = tmp_path / "pipeline_log.md"
        log.write_text("x" * (3 * _LOG_LINE_MAX_CHARS + 7) + "\nnext\n")
        manager.get(job_id).log_file = str(log)

        body = client.get(f"/v1/jobs/{job_id}/log").json()
        assert len(body["lines"][0]) == _LOG_LINE_MAX_CHARS + len(_LOG_LINE_TRUNCATED_SUFFIX)
        assert body["lines"][0].endswith(_LOG_LINE_TRUNCATED_SUFFIX)
        assert body["lines"][1] == "next"


class TestJobLogWindow:
    def test_reads_only_the_requested_window(self):
        import io

        from refgenie.server.jobs.router import _log_window

        content = "".join(f"line {i}\n" for i in range(200_000))
        handle = io.StringIO(content)

        window, has_more = _log_window(handle, offset=10, limit=5)

        assert window == [f"line {i}" for i in range(10, 15)]
        assert has_more is True
        assert handle.tell() < len(content) * 0.01

    def test_first_page(self):
        import io

        from refgenie.server.jobs.router import _log_window

        content = "".join(f"line {i}\n" for i in range(200_000))
        handle = io.StringIO(content)

        window, has_more = _log_window(handle, offset=0, limit=5)

        assert window == [f"line {i}" for i in range(5)]
        assert has_more is True
        assert handle.tell() < len(content) * 0.01


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

    def test_delete_forgets_a_terminal_job_and_refuses_a_running_one(
        self, job_manager_factory, gate
    ):
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
        for path in (
            "/v1/jobs",
            f"/v1/jobs/{job_id}",
            f"/v1/jobs/{job_id}/events",
            f"/v1/jobs/{job_id}/log",
        ):
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
            until=lambda f: (
                f["event"] == "done"
                and later.get("ref") is not None
                and f["data"]["job_id"] == later["ref"].job_id
            ),
            during=submit_a_second_job,
        )
        done_ids = [f["data"]["job_id"] for f in frames if f["event"] == "done"]
        assert done_ids == [first.job_id, later["ref"].job_id]

    def test_last_event_id_and_since_resume_from_the_global_sequence(self, job_manager_factory):
        manager = job_manager_factory()
        job = manager.submit_pull(pull_params())
        manager.wait(job.job_id, timeout=5)

        _, _, by_header = read_sse(manager, until=until_done, headers=[(b"last-event-id", b"2")])
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
            "id",
            "kind",
            "status",
            "label",
            "target",
            "queue_position",
            "progress",
            "created_at",
            "started_at",
            "finished_at",
            "result",
            "error",
            "log_lines",
            "cancellable",
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
            "status",
            "progress",
            "log",
            "done",
            "truncated",
            "heartbeat",
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
