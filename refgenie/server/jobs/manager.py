"""The in-process job manager: runs pulls and builds on dedicated thread pools.

Jobs live in memory only and do not survive a restart. See "Background jobs in
local mode" in docs/design-notes.md.
"""

import hashlib
import json
import logging
import threading
import time
import traceback
import uuid
from collections import deque
from collections.abc import Callable, Collection, Mapping
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime
from enum import StrEnum
from typing import TYPE_CHECKING

from refgenie import progress
from refgenie.logger import logger
from refgenie.progress import ProgressAborted, ProgressEvent
from refgenie.server.errors import ErrorCode, classify_exception
from refgenie.server.jobs.logcapture import JobLogHandler
from refgenie.server.jobs.logcapture import install as install_log_handler
from refgenie.server.jobs.logcapture import remove as remove_log_handler
from refgenie.server.jobs.schemas import (
    TERMINAL_STATUSES,
    BuildJobParams,
    GenomeInitJobParams,
    JobError,
    JobEvent,
    JobKind,
    JobLinks,
    JobParams,
    JobProgress,
    JobRecord,
    JobRef,
    JobResult,
    JobStatus,
    JobTarget,
    PullJobParams,
    describe,
    target_of,
    utcnow,
)

if TYPE_CHECKING:
    from refgenie.core import Refgenie

__all__ = ["CancelOutcome", "JobContext", "JobManager", "Runner"]


class CancelOutcome(StrEnum):
    """The three answers `request_cancel` can give."""

    ACCEPTED = "accepted"
    NOT_CANCELLABLE = "not_cancellable"
    ALREADY_DONE = "already_done"


@dataclass
class JobContext:
    """What a runner is handed: reporting and cancellation hooks.

    Runners never set status; the manager owns the state machine.
    """

    job_id: str
    params: JobParams
    refgenie: "Refgenie"
    #: Announce a coarse phase change: phase("download", "Downloading rCRSd").
    #: Names come from `frontend/src/components/jobs/phases.ts`.
    phase: Callable[..., None]
    #: Append one log line to this job.
    log: Callable[..., None]
    #: Raises ProgressAborted if cancellation was requested. Call it in any
    #: loop that does not already go through a progress sink.
    check_cancel: Callable[[], None]
    #: Record the pypiper log path for a build, so the UI can link to it.
    set_log_file: Callable[[str], None]


Runner = Callable[[JobContext], "JobResult | None"]

#: How many pulls may run at once. Safe because no `rich` live display is built
#: while a progress sink is installed.
DEFAULT_PULL_WORKERS = 2


# ---------------------------------------------------------------------------
# Internal records (not on the wire; guarded by the manager's lock)
# ---------------------------------------------------------------------------


@dataclass
class Job:
    """One background job, mutated in place under the manager's lock.

    Readers get a `JobRecord` snapshot instead.
    """

    id: str
    kind: JobKind
    params: JobParams
    #: Stable hash of the parameters, for duplicate coalescing.
    params_key: str
    label: str
    target: JobTarget
    created_at: datetime = field(default_factory=utcnow)
    status: JobStatus = JobStatus.QUEUED
    started_at: datetime | None = None
    finished_at: datetime | None = None
    progress: JobProgress | None = None
    result: JobResult | None = None
    error: JobError | None = None
    #: Count of `log` events appended for this job.
    log_lines: int = 0
    skipped: bool | None = None
    cancel_requested: threading.Event = field(default_factory=threading.Event)
    #: pypiper's log path for builds; the UI links to it and tails it.
    log_file: str | None = None
    #: Highest manager-global seq assigned to an event of this job.
    event_seq: int = 0

    @property
    def is_terminal(self) -> bool:
        return self.status in TERMINAL_STATUSES


@dataclass
class EventRing:
    """The event log for all jobs: a bounded deque plus a global seq.

    `dropped_before` tracks overflow so a resuming client learns what it missed.
    """

    maxlen: int
    events: deque[JobEvent] = field(init=False)
    last_seq: int = 0
    dropped_before: int = 1

    def __post_init__(self) -> None:
        self.events = deque(maxlen=self.maxlen)

    def append(self, make_event) -> JobEvent:
        """Assign the next seq, build the event and store it.

        Args:
            make_event: Callable taking the assigned seq and returning a
                `JobEvent`.
        """
        self.last_seq += 1
        event = make_event(self.last_seq)
        if len(self.events) == self.maxlen and self.events:
            self.dropped_before = self.events[0].seq + 1
        self.events.append(event)
        return event

    def since(self, seq: int) -> tuple[list[JobEvent], int, int]:
        """Events after `seq`, the new cursor, and how many unseen events were dropped."""
        dropped = max(0, self.dropped_before - (seq + 1))
        selected = [event for event in self.events if event.seq > seq]
        cursor = selected[-1].seq if selected else max(seq, self.last_seq)
        return selected, cursor, dropped


class JobManager:
    """Owns every background job in this process.

    Args:
        refgenie: The `Refgenie` instance runners operate on.
        runners: Override the runner for a kind. Tests pass fakes here; in
            production the defaults from `refgenie.server.jobs.runners` are used.
        pull_workers: Concurrent pull slots.
        history_limit: How many jobs to retain. Only finished jobs are evicted.
        event_buffer: How many events the manager-global ring holds.
    """

    def __init__(
        self,
        refgenie: "Refgenie",
        *,
        runners: Mapping[JobKind, Runner] | None = None,
        pull_workers: int = DEFAULT_PULL_WORKERS,
        history_limit: int = 100,
        event_buffer: int = 2000,
    ) -> None:
        self._refgenie = refgenie
        self._history_limit = history_limit

        if runners is None:
            from refgenie.server.jobs import runners as default_runners

            runners = default_runners.DEFAULT_RUNNERS
        self._runners: dict[JobKind, Runner] = dict(runners)

        self._lock = threading.RLock()
        self._jobs: dict[str, Job] = {}
        self._order: list[str] = []
        self._futures: dict[str, Future] = {}
        self._ring = EventRing(maxlen=event_buffer)
        #: thread ident -> job id, for attributing log records.
        self._thread_jobs: dict[int, str] = {}
        #: Per-executor FIFO of job ids not yet started, for queue_position.
        self._queues: dict[str, list[str]] = {"pull": [], "build": []}
        self._shutdown = False

        self._pull_pool = ThreadPoolExecutor(
            max_workers=pull_workers, thread_name_prefix="refgenie-pull"
        )
        # One build at a time: pypiper swaps sys.stdout/sys.stderr process-wide,
        # and overlapping builds corrupt that save/restore chain.
        self._build_pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="refgenie-build")

        # Each manager's handler only recognizes its own worker threads.
        self._log_handler = JobLogHandler(self._route_log, level=logging.INFO)
        install_log_handler(self._log_handler)

    @property
    def refgenie(self) -> "Refgenie":
        """The instance every runner operates on."""
        return self._refgenie

    # -- submission ---------------------------------------------------------

    def submit_pull(self, params: PullJobParams) -> JobRef:
        """Queue a pull. Never blocks; never raises for a job-level failure."""
        return self._submit(JobKind.PULL, params)

    def submit_build(self, params: BuildJobParams) -> JobRef:
        """Queue a build. Never blocks; never raises for a job-level failure."""
        return self._submit(JobKind.BUILD, params)

    def submit_genome_init(self, params: GenomeInitJobParams) -> JobRef:
        """Queue a genome initialization (and its fasta build)."""
        return self._submit(JobKind.GENOME_INIT, params)

    def _submit(self, kind: JobKind, params: JobParams) -> JobRef:
        params_key = _params_key(kind, params)
        with self._lock:
            if self._shutdown:
                raise RuntimeError("JobManager is shut down; no new jobs accepted.")
            # Coalesce a duplicate onto the active job instead of starting a second.
            for existing in self._jobs.values():
                if existing.params_key == params_key and not existing.is_terminal:
                    return self._ref(existing, duplicate=True)

            job = Job(
                id=uuid.uuid4().hex[:12],
                kind=kind,
                params=params,
                params_key=params_key,
                label=describe(kind, params),
                target=target_of(kind, params),
            )
            self._jobs[job.id] = job
            self._order.append(job.id)
            self._queues[_pool_key(kind)].append(job.id)
            self._append_event(
                job,
                "status",
                status=JobStatus.QUEUED,
                message="Queued",
                queue_position=self._queue_position(job),
            )
            ref = self._ref(job)

        pool = self._build_pool if _pool_key(kind) == "build" else self._pull_pool
        future = pool.submit(self._run, job.id)
        with self._lock:
            self._futures[job.id] = future
        return ref

    # -- reads --------------------------------------------------------------

    def record(self, job_id: str) -> JobRecord:
        """An immutable snapshot of one job. Raises KeyError if unknown."""
        with self._lock:
            return self._record(self._jobs[job_id])

    def list(
        self,
        *,
        kind: JobKind | None = None,
        status: "JobStatus | Collection[JobStatus] | None" = None,
    ) -> list[JobRecord]:
        """Job snapshots, newest first, filtered by kind and one or more statuses."""
        with self._lock:
            records = [
                self._record(self._jobs[job_id])
                for job_id in reversed(self._order)
                if job_id in self._jobs
            ]
        if kind is not None:
            records = [r for r in records if r.kind == kind]
        if status is not None:
            wanted = {status} if isinstance(status, JobStatus) else set(status)
            records = [r for r in records if r.status in wanted]
        return records

    # Quoted: the `list` method shadows the builtin in this class body.
    def events_since(self, seq: int) -> "tuple[list[JobEvent], int, int]":
        """Events after `seq` across all jobs, the new cursor, and how many unseen events were dropped."""
        with self._lock:
            return self._ring.since(seq)

    def wait(self, job_id: str, timeout: float = 10.0) -> JobRecord:
        """Block until the job reaches a terminal state, then return it.

        Polls status rather than the future, since a job cancelled while queued
        has no future to wait on.

        Raises:
            TimeoutError: If the job is still running when `timeout` expires.
        """
        deadline = time.monotonic() + timeout
        while True:
            with self._lock:
                job = self._jobs[job_id]
                if job.is_terminal:
                    return self._record(job)
            if time.monotonic() >= deadline:
                raise TimeoutError(f"Job {job_id} did not finish within {timeout}s")
            time.sleep(0.01)

    # -- cancellation and cleanup ------------------------------------------

    def request_cancel(self, job_id: str) -> CancelOutcome:
        """Ask a job to stop. Raises KeyError if unknown.

        A queued job is cancelled outright and a running pull stops at its next
        progress report. A running build cannot be cancelled: pypiper exposes no
        handle to its subprocesses.
        """
        with self._lock:
            job = self._jobs[job_id]
            if job.is_terminal:
                return CancelOutcome.ALREADY_DONE
            job.cancel_requested.set()
            if job.status == JobStatus.QUEUED:
                future = self._futures.get(job_id)
                if future is not None and future.cancel():
                    self._queue_remove(job)
                    self._finish(job, JobStatus.CANCELLED, error=_cancelled_error())
                # Otherwise it is already dispatched; `_run` re-checks the flag.
                return CancelOutcome.ACCEPTED
            if job.kind in (JobKind.BUILD, JobKind.GENOME_INIT):
                return CancelOutcome.NOT_CANCELLABLE
            return CancelOutcome.ACCEPTED

    def forget(self, job_id: str) -> None:
        """Drop a terminal job from the registry. Raises KeyError if unknown,
        ValueError if it is still running."""
        with self._lock:
            job = self._jobs[job_id]
            if not job.is_terminal:
                raise ValueError(f"Job {job_id} is not finished.")
            self._drop(job_id)

    def shutdown(self, wait: bool = False) -> None:
        """Stop accepting work. By default does not wait; abandoned jobs are logged."""
        with self._lock:
            if self._shutdown:
                return
            self._shutdown = True
            running = [
                job.id
                for job in self._jobs.values()
                if job.status in (JobStatus.RUNNING, JobStatus.QUEUED)
            ]
        if running:
            logger.warning(
                f"{len(running)} job(s) abandoned at process exit ({', '.join(running)}); "
                "jobs are in-memory and do not survive a restart."
            )
        self._pull_pool.shutdown(wait=wait, cancel_futures=True)
        self._build_pool.shutdown(wait=wait, cancel_futures=True)
        remove_log_handler(self._log_handler)

    # -- the executor callable ---------------------------------------------

    def _run(self, job_id: str) -> None:
        """Run one job. The only place that moves a job out of `queued`."""
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None or job.is_terminal:
                return
            if job.cancel_requested.is_set():
                self._queue_remove(job)
                self._finish(job, JobStatus.CANCELLED, error=_cancelled_error())
                return
            ident = threading.get_ident()
            job.status = JobStatus.RUNNING
            job.started_at = utcnow()
            self._queue_remove(job)
            self._thread_jobs[ident] = job.id
            self._append_event(job, "status", status=JobStatus.RUNNING, message="Running")
            runner = self._runners[job.kind]
            ctx = JobContext(
                job_id=job.id,
                params=job.params,
                refgenie=self._refgenie,
                phase=lambda name, message=None: self._set_phase(job.id, name, message),
                log=lambda line, source="refgenie": self._emit_from_runner(
                    job.id, "log", line=line, source=source
                ),
                check_cancel=lambda: self._check_cancel(job.id),
                set_log_file=lambda path: self._set_log_file(job.id, path),
            )

        status = JobStatus.SUCCEEDED
        result: JobResult | None = None
        error: JobError | None = None
        try:
            with progress.use_sink(self._make_sink(job.id)):
                result = runner(ctx)
        except ProgressAborted:
            status = JobStatus.CANCELLED
            error = _cancelled_error()
        except Exception as exc:  # noqa: BLE001 - the manager IS the boundary
            _, code, message = classify_exception(exc)
            # A skipped pull declined to overwrite; report it as cancelled, not failed.
            status = JobStatus.CANCELLED if code == ErrorCode.PULL_SKIPPED else JobStatus.FAILED
            error = JobError(code=code, message=message, detail=traceback.format_exc())
            logger.debug(f"Job {job.id} ({job.kind}) {status}: {message}")
        finally:
            with self._lock:
                self._thread_jobs.pop(threading.get_ident(), None)
                self._finish(job, status, result=result, error=error)
                self._trim_history()

    # -- internals ----------------------------------------------------------

    def _finish(
        self,
        job: Job,
        status: JobStatus,
        *,
        result: JobResult | None = None,
        error: JobError | None = None,
    ) -> None:
        """Move a job to a terminal state and announce it. Caller holds the lock."""
        if job.is_terminal:
            return
        job.status = status
        job.result = result
        job.error = error
        job.finished_at = utcnow()
        if job.kind in (JobKind.BUILD, JobKind.GENOME_INIT):
            # A succeeded build with no pypiper log never ran: the asset existed.
            job.skipped = status == JobStatus.SUCCEEDED and job.log_file is None
        self._append_event(job, "status", status=status, message=status.value)
        self._append_event(job, "done", status=status)

    def _drop(self, job_id: str) -> None:
        """Remove a job from the registry. Caller holds the lock."""
        self._jobs.pop(job_id, None)
        self._futures.pop(job_id, None)
        if job_id in self._order:
            self._order.remove(job_id)

    def _trim_history(self) -> None:
        """Evict the oldest finished jobs past the limit. Caller holds the lock."""
        while len(self._order) > self._history_limit:
            for job_id in list(self._order):
                job = self._jobs.get(job_id)
                if job is None or job.is_terminal:
                    self._drop(job_id)
                    break
            else:
                return

    def _queue_remove(self, job: Job) -> None:
        queue = self._queues[_pool_key(job.kind)]
        if job.id in queue:
            queue.remove(job.id)

    def _queue_position(self, job: Job) -> int | None:
        if job.status != JobStatus.QUEUED:
            return None
        queue = self._queues[_pool_key(job.kind)]
        return queue.index(job.id) if job.id in queue else None

    def _ref(self, job: Job, duplicate: bool = False) -> JobRef:
        return JobRef(
            job_id=job.id,
            kind=job.kind,
            status=job.status,
            created_at=job.created_at,
            duplicate=duplicate,
            links=self._links(job.id),
        )

    def _links(self, job_id: str) -> JobLinks:
        return JobLinks(
            self=f"/v1/jobs/{job_id}",
            events="/v1/jobs/events",
            cancel=f"/v1/jobs/{job_id}/cancel",
        )

    def _record(self, job: Job) -> JobRecord:
        """Snapshot a job. Caller holds the lock."""
        return JobRecord(
            id=job.id,
            kind=job.kind,
            status=job.status,
            label=job.label,
            target=job.target,
            created_at=job.created_at,
            links=self._links(job.id),
            params=job.params.model_dump(mode="json"),
            started_at=job.started_at,
            finished_at=job.finished_at,
            queue_position=self._queue_position(job),
            progress=job.progress,
            result=job.result,
            error=job.error,
            cancellable=_cancellable(job),
            log_lines=job.log_lines,
            skipped=job.skipped,
            event_seq=job.event_seq,
            log_file=job.log_file,
        )

    def _append_event(self, job: Job, type: str, **fields) -> JobEvent:
        """Append one event to the manager-global ring. Caller holds the lock."""

        def make(seq: int) -> JobEvent:
            return JobEvent(job_id=job.id, seq=seq, ts=utcnow(), type=type, **fields)

        event = self._ring.append(make)
        job.event_seq = event.seq
        if type == "log":
            job.log_lines += 1
        return event

    def _emit_from_runner(self, job_id: str, type: str, **fields) -> None:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                return
            self._append_event(job, type, **fields)

    def _set_phase(self, job_id: str, name: str, message: str | None = None) -> None:
        """Record a coarse phase change and announce it as a progress frame."""
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                return
            job.progress = JobProgress(phase=name, message=message)
            self._append_event(job, "progress", phase=name, message=message)

    def _set_log_file(self, job_id: str, path: str) -> None:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is not None:
                job.log_file = path

    def _check_cancel(self, job_id: str) -> None:
        with self._lock:
            job = self._jobs.get(job_id)
            requested = job is not None and job.cancel_requested.is_set()
        if requested:
            raise ProgressAborted("Cancelled by user.")

    def _make_sink(self, job_id: str):
        """The progress sink installed for one job's whole run.

        Raising `ProgressAborted` here is how a running pull is cancelled; the
        pull transaction rolls back on the way out.
        """

        def sink(event: ProgressEvent) -> None:
            with self._lock:
                job = self._jobs.get(job_id)
                if job is None:
                    return
                if job.cancel_requested.is_set():
                    raise ProgressAborted("Cancelled by user.")

                if event.type == "log":
                    self._append_event(job, "log", line=event.message, source="refgenie")
                    return

                # "stage" and "progress" library events both become `progress` frames.
                job.progress = _merge_progress(job.progress, event)
                self._append_event(
                    job,
                    "progress",
                    phase=job.progress.phase,
                    percent=job.progress.percent,
                    bytes_done=job.progress.bytes_done,
                    bytes_total=job.progress.bytes_total,
                    message=job.progress.message,
                )

        return sink

    def _route_log(self, thread_id: int, level: str, message: str) -> None:  # noqa: D401
        """Attribute one `refgenie` log record to the job on `thread_id`."""
        with self._lock:
            job_id = self._thread_jobs.get(thread_id)
            if job_id is None:
                return
            job = self._jobs.get(job_id)
            if job is None:
                return
            self._append_event(job, "log", line=message, level=level, source="refgenie")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _pool_key(kind: JobKind) -> str:
    """Which executor runs `kind`; genome_init ends in a build, so it queues with builds."""
    return "pull" if kind == JobKind.PULL else "build"


def _cancellable(job: Job) -> bool:
    """Whether cancelling would do anything: queued jobs and running pulls only."""
    if job.is_terminal:
        return False
    if job.status == JobStatus.QUEUED:
        return True
    return job.kind == JobKind.PULL


def _cancelled_error() -> JobError:
    return JobError(code=ErrorCode.CANCELLED, message="Cancelled by user.")


def _percent(current: int | None, total: int | None) -> float | None:
    """0-100, or None when the total is unknown (an indeterminate bar)."""
    if current is None or not total:
        return None
    return min(100.0, round(current / total * 100, 1))


def _merge_progress(existing: JobProgress | None, event: ProgressEvent) -> JobProgress:
    """Fold one library progress event into the job's progress reading.

    Only the phase carries over, when the event has none; percent, byte counts
    and message come from the event alone and are None when it lacks them.
    """
    phase = event.extra.get("phase") or (existing.phase if existing else None)
    bytes_done = bytes_total = None
    if event.unit == "bytes":
        bytes_done, bytes_total = event.current, event.total
    return JobProgress(
        phase=phase or "working",
        percent=_percent(event.current, event.total),
        message=event.message,
        bytes_done=bytes_done,
        bytes_total=bytes_total,
    )


def _params_key(kind: JobKind, params: JobParams) -> str:
    """A stable hash of (kind, params), for duplicate coalescing."""
    payload = json.dumps(
        {"kind": str(kind), "params": params.model_dump(mode="json")},
        sort_keys=True,
        default=str,
    )
    return hashlib.sha256(payload.encode()).hexdigest()
