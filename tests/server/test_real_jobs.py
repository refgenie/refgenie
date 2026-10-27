"""Real pulls and builds through the job machinery.

Where ``test_jobs.py`` drives ``JobManager`` with fake runners, these prove
the seams meet with real work: a build from a worker thread, stage validation
before any build, one real pull all the way through the manager, and the
assembled ``refgenie dash`` app from boot to a terminal job over HTTP.

Tier: the classes that build or pull real assets carry a class-level
``pytest.mark.component``; the stage-validation class uses a test double and
says ``unit`` explicitly.
"""

import threading
from pathlib import Path

import pytest

from tests.helpers import (
    ACTION_HEADERS,
    ASSET,
    GENOME,
    GROUP,
    PHASE_VOCABULARY,
    make_local_client,
    requires_dash,
    requires_server,
)

requires_server()
requires_dash()

from refgenie.server.errors import ErrorCode  # noqa: E402  (must follow the extras guard)
from refgenie.server.jobs.manager import JobManager  # noqa: E402
from refgenie.server.jobs.schemas import JobStatus, PullJobParams  # noqa: E402


# ===========================================================================
# Building from a worker thread; failing fast on a bad stage request
# ===========================================================================


class TestBuildOnWorkerThread:
    """The build registers its SIGINT handler only on the main thread:
    ``signal.signal()`` raises off the main thread, so an unconditional handler
    kills every threaded build on its first line of real work. Both tests build
    real assets: component."""

    pytestmark = pytest.mark.component

    def test_build_succeeds_on_a_worker_thread(self, refgenie_fs):
        """The anchor test for the main-thread guard on ``signal.signal``."""
        box = {}

        def build():
            try:
                box["asset"] = refgenie_fs.build.run(
                    recipe_name="fasta",
                    genome_alias=GENOME,
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
        import refgenie.managers.build as builder_module

        registered = []
        real_signal = builder_module.signal.signal

        def spy(sig, handler):
            registered.append(sig)
            return real_signal(sig, handler)

        monkeypatch.setattr(builder_module.signal, "signal", spy)
        refgenie_fs.build.run(
            recipe_name="fasta", genome_alias=GENOME, asset_group_name=GROUP, asset_name=ASSET
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
        monkeypatch.setattr(unstaged.build, "_build", lambda *a, **kw: called.append(kw))

        with pytest.raises(ValueError, match="genome_stage_folder is not set"):
            unstaged.build.run(
                recipe_name="fasta",
                genome_alias=GENOME,
                asset_group_name=GROUP,
                asset_name=ASSET,
                stage=True,
            )
        assert called == [], "the build ran before the stage folder was validated"

    def test_no_staging_requested_is_unaffected(self, unstaged, monkeypatch):
        monkeypatch.setattr(unstaged.alias, "resolve", lambda name: "d" * 64)
        monkeypatch.setattr(unstaged.build, "_build", lambda *a, **kw: None)

        assert (
            unstaged.build.run(
                recipe_name="fasta", genome_alias=GENOME, asset_group_name=GROUP, asset_name=ASSET
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
