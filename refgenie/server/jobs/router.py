"""The jobs HTTP API, including the one multiplexed SSE stream.

Mounted at ``/v1`` in local mode, so the paths below read ``/v1/jobs/...``.
Nothing here imports ``create_app``: the router takes its manager off
``request.app.state``, so a test can include it into a bare ``FastAPI()`` and
exercise the real routes.

SSE payload contract -- the frontend codes against this verbatim
-----------------------------------------------------------------

There is exactly ONE stream, ``GET /v1/jobs/events``, multiplexed across every
job. Not one per job: browsers cap HTTP/1.1 at roughly six connections per
origin, and ``EventSource`` cannot set request headers, so a stream-per-job is
both wasteful and unauthenticable.

::

    id: 17
    event: progress
    data: {"job_id":"a1b2c3d4e5f6","seq":17,"ts":"2026-08-12T18:22:03.412Z",
           "type":"progress","phase":"download","percent":30.0,
           "bytes_done":12582912,"bytes_total":41943040,
           "message":"rCRSd/fasta"}

    id: 18
    event: log
    data: {"job_id":"...","seq":18,"ts":"...","type":"log","level":"INFO",
           "source":"refgenie","line":"Extracting asset tarball: ..."}

    id: 42
    event: done
    data: {<the full JobRecord for ONE job, keyed `id`>}

    id: 43
    event: heartbeat
    data: {"seq":43,"ts":"..."}

Rules a client may rely on:

* ``event:`` is one of exactly six names -- ``status``, ``progress``, ``log``,
  ``done``, ``truncated``, ``heartbeat``. There is no seventh, and adding one
  would not be a compatible extension: ``EventSource`` dispatches by name and
  the client registers a listener per known name, so an unrecognized event is
  dropped by the browser before any JS sees it. A coarse phase change is a
  ``progress`` frame with a ``phase`` and no byte counts.
* ``data:`` is always a single-line JSON object, serialized without its null
  fields -- a ``log`` frame does not carry four empty progress keys.
* Every event carries a ``job_id`` except ``truncated``, which describes the
  manager's ring rather than any one job.
* A ``log`` frame carries ``line`` (not ``message``) and a ``source`` of
  ``refgenie`` or ``pipeline``.
* A ``progress`` frame carries the ``JobProgress`` fields flattened to the top
  level: ``phase``, ``percent``, ``message``, ``bytes_done``, ``bytes_total``.
  The client merges them into the job's cumulative reading.
* ``id:`` is the event ``seq``, **manager-global and monotonic across all
  jobs**. Reconnection sends ``Last-Event-ID``; ``?since=<seq>`` is the
  explicit form of the same thing.
* Exactly one ``done`` per job, carrying that job's complete ``JobRecord``.
  The stream stays open afterwards. There are no separate ``result`` or
  ``error`` events -- branch on ``JobRecord.status``.
* The keepalive is a named ``heartbeat`` event every 15 s, never a
  ``: keepalive`` comment: ``EventSource`` does not dispatch comments to JS,
  so a comment would leave a client watchdog with nothing to observe. The
  client's watchdog fires at three missed beats, so this interval is a promise.
* A ``truncated`` frame means the ring dropped events the client had not seen;
  it carries ``dropped`` (how many). Re-fetch ``GET /v1/jobs`` for
  authoritative state.
* No custom header is required on this route, and none is accepted --
  ``EventSource`` cannot send one.

Jobs are in-process memory. Restarting `refgenie dash` loses every job record;
a pull interrupted that way may leave a partial ``.tgz`` in the asset group
directory, which the next pull overwrites.
"""

import asyncio
import json
from typing import Any

from fastapi import APIRouter, Depends, Header, HTTPException, Query, Request, Response
from fastapi.responses import StreamingResponse

from refgenie.const import DEFAULT_PAGE_SIZE, MAX_PAGE_SIZE
from refgenie.server.errors import ErrorCode
from refgenie.server.jobs.manager import CancelOutcome, JobManager
from refgenie.server.jobs.schemas import (
    PARAMS_MODEL,
    BuildJobParams,
    JobEvent,
    JobKind,
    JobRecord,
    JobRef,
    JobStatus,
    JobStatusFilter,
    JobSubmission,
    utcnow,
)
from refgenie.utils.pagination import PaginatedResponse, PaginationMeta

__all__ = ["get_job_manager", "router"]

router = APIRouter(prefix="/jobs", tags=["Jobs"])


def _require_action_header(request: Request) -> None:
    """The same anti-CSRF guard the actions router carries, lazily imported.

    Lazy on purpose: ``refgenie.server.local`` imports this package (its
    ``__init__`` pulls in the actions router, which imports ``get_job_manager``
    from here), so a module-level import of ``local.security`` would be a
    circular import. The guard itself is identical -- one owner, in
    ``refgenie/server/local/security.py``.
    """
    from refgenie.server.local.security import require_action_header

    require_action_header(request)


def _require_action_origin(request: Request) -> None:
    """See :func:`_require_action_header`; the bridge origin policy layer."""
    from refgenie.server.local.security import require_action_origin

    require_action_origin(request)


#: Attached to the three state-changing routes (submit, cancel, forget) --
#: NEVER to the reads, and never to the SSE stream (``EventSource`` cannot set
#: request headers). The read surface is deliberately open to bridge origins.
_WRITE_GUARDS = [Depends(_require_action_header), Depends(_require_action_origin)]

#: How often the SSE generator looks for new events. 250 ms is imperceptible
#: against a download that takes minutes, and polling an in-memory ring from an
#: `async def` costs nothing -- unlike pushing from worker threads through an
#: asyncio.Queue, which would need a cross-thread loop handle and would still
#: leave the JSON polling fallback to be written separately.
POLL_SECONDS = 0.25
#: Seconds of quiet before a heartbeat frame.
HEARTBEAT_SECONDS = 15.0


def get_job_manager(request: Request) -> JobManager:
    """The app's JobManager, or a 503 in server mode.

    Server mode runs many workers behind a load balancer and manages nobody's
    assets, so it has no manager to hand out. 503 rather than 404 says "not
    here", not "you got the URL wrong".
    """
    manager = getattr(request.app.state, "job_manager", None)
    if manager is None:
        raise HTTPException(
            status_code=503, detail="Background jobs are available in local mode only."
        )
    return manager


def _status_filter(status: "JobStatusFilter | None") -> "frozenset[JobStatus] | None":
    """Expand a status query parameter into the concrete statuses it selects."""
    return None if status is None else status.statuses()


def _job_or_404(manager: JobManager, job_id: str) -> JobRecord:
    try:
        return manager.record(job_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"No such job: {job_id}") from None


# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------


@router.post("", status_code=201, response_model=JobRef, dependencies=_WRITE_GUARDS)
def submit_job(
    submission: JobSubmission,
    response: Response,
    manager: JobManager = Depends(get_job_manager),
) -> JobRef:
    """Submit a job of any kind.

    The actions router has its own typed endpoints and calls
    `manager.submit_*` directly; this generic form exists so the jobs API is
    complete and drivable by hand.
    """
    model = PARAMS_MODEL[submission.kind]
    try:
        params = model(**submission.params)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    _preflight(manager, submission.kind, params)

    if submission.kind == JobKind.PULL:
        ref = manager.submit_pull(params)
    elif submission.kind == JobKind.BUILD:
        ref = manager.submit_build(params)
    else:
        ref = manager.submit_genome_init(params)
    response.headers["Location"] = ref.links.self
    return ref


def _preflight(manager: JobManager, kind: JobKind, params) -> None:
    """Checks worth doing synchronously, so the user learns in milliseconds.

    A cheap backstop only. `Refgenie.build_asset` validates the same thing
    first thing it does; this exists so a misconfigured stage folder is a 400
    on the click rather than a failed job card.
    """
    if kind == JobKind.BUILD and isinstance(params, BuildJobParams) and params.stage:
        if getattr(manager.refgenie, "genome_stage_folder", None) is None:
            raise HTTPException(
                status_code=400,
                detail=(
                    "Staging is not configured: genome_stage_folder is not set. "
                    "Set it in the refgenie config or build without staging."
                ),
            )


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------


@router.get("", response_model=PaginatedResponse[JobRecord])
def list_jobs(
    kind: JobKind | None = None,
    status: JobStatusFilter | None = None,
    offset: int = Query(0, ge=0),
    limit: int = Query(DEFAULT_PAGE_SIZE, ge=1, le=MAX_PAGE_SIZE),
    manager: JobManager = Depends(get_job_manager),
) -> PaginatedResponse[JobRecord]:
    """Job records, newest first, in the same envelope as every other listing.

    `status` also takes two pseudo-values, because they are what a console
    actually asks for: `active` is queued+running and `terminal` is
    succeeded/failed/cancelled. The UI polls for both on every cycle, and
    spelling that as three separate requests would triple the traffic.
    """
    records = manager.list(kind=kind, status=_status_filter(status))
    window = records[offset : offset + limit]
    return PaginatedResponse[JobRecord](
        items=window,
        pagination=PaginationMeta(offset=offset, limit=limit, total=len(records)),
    )


# Registered BEFORE /{job_id} so "events" is never captured as a job id.
@router.get("/events")
async def stream_events(
    request: Request,
    since: int = 0,
    last_event_id: str | None = Header(None, alias="Last-Event-ID"),
    manager: JobManager = Depends(get_job_manager),
) -> StreamingResponse:
    """The one multiplexed SSE stream. See the module docstring for the contract.

    A stream opened *after* a job finished still replays from the cursor and
    still delivers that job's `done`, which removes the submit-then-connect
    race: a client may connect whenever it likes.
    """
    cursor = _resume_cursor(since, last_event_id)

    async def generate():
        nonlocal cursor
        quiet = 0.0
        while True:
            if await request.is_disconnected():
                return
            events, cursor, dropped = manager.events_since_detailed(cursor)
            if dropped:
                yield _frame(
                    "truncated",
                    {
                        "seq": cursor,
                        "ts": utcnow().isoformat(),
                        "type": "truncated",
                        "dropped": dropped,
                        "message": "Older events were dropped; re-fetch /v1/jobs.",
                    },
                )
            for event in events:
                yield _frame(event.type, _event_payload(manager, event), seq=event.seq)
            if events:
                quiet = 0.0
                continue
            await asyncio.sleep(POLL_SECONDS)
            quiet += POLL_SECONDS
            if quiet >= HEARTBEAT_SECONDS:
                quiet = 0.0
                yield _frame(
                    "heartbeat", {"seq": cursor, "ts": utcnow().isoformat()}, seq=cursor
                )

    return StreamingResponse(
        generate(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            # nginx buffers proxied responses by default, which would hold a
            # progress stream until the download finished.
            "X-Accel-Buffering": "no",
        },
    )


@router.get("/{job_id}", response_model=JobRecord)
def get_job(job_id: str, manager: JobManager = Depends(get_job_manager)) -> JobRecord:
    """One job record."""
    return _job_or_404(manager, job_id)


@router.get("/{job_id}/events", response_model=dict)
def poll_job_events(
    job_id: str,
    since: int = 0,
    manager: JobManager = Depends(get_job_manager),
) -> dict[str, Any]:
    """Polling fallback for one job. There is no per-job SSE.

    Returns this job's slice of the manager-global log, plus the cursor to
    resume from and the job's current record.
    """
    record = _job_or_404(manager, job_id)
    events, cursor, truncated = manager.events_since(since)
    return {
        # Serialized by the same function the stream uses, so a client that
        # falls back to polling parses identical frames -- including a `done`
        # that carries the whole record.
        "events": [
            _event_payload(manager, event) for event in events if event.job_id == job_id
        ],
        "next": cursor,
        "truncated": truncated,
        "job": record.model_dump(mode="json"),
    }


@router.get("/{job_id}/log", response_model=dict)
def get_job_log(
    job_id: str,
    offset: int = 0,
    limit: int = Query(500, ge=1, le=5000),
    manager: JobManager = Depends(get_job_manager),
) -> dict[str, Any]:
    """Lines from the job's pipeline log file, from `offset`.

    `offset` is a line index, not a byte offset, so a client can page forward
    without knowing anything about encoding.
    """
    record = _job_or_404(manager, job_id)
    if not record.log_file:
        return {"lines": [], "next_offset": offset, "truncated": False}
    try:
        with open(record.log_file, "r", errors="replace") as handle:
            all_lines = handle.read().splitlines()
    except OSError:
        return {"lines": [], "next_offset": offset, "truncated": False}
    window = all_lines[offset : offset + limit]
    return {
        "lines": window,
        "next_offset": offset + len(window),
        "truncated": offset + len(window) < len(all_lines),
    }


# ---------------------------------------------------------------------------
# Writes
# ---------------------------------------------------------------------------


@router.post("/{job_id}/cancel", status_code=202, response_model=JobRecord, dependencies=_WRITE_GUARDS)
def cancel_job(job_id: str, manager: JobManager = Depends(get_job_manager)) -> JobRecord:
    """Ask a job to stop.

    409 with a distinguishing code when it cannot: a running build owns a
    subprocess tree pypiper will not give up, and a finished job has nothing
    left to stop.
    """
    _job_or_404(manager, job_id)
    outcome = manager.request_cancel(job_id)
    if outcome == CancelOutcome.NOT_CANCELLABLE:
        raise HTTPException(
            status_code=409,
            detail={
                "code": "not_cancellable",
                "message": (
                    "A running build cannot be cancelled. Stop the refgenie dash "
                    "process to abort it."
                ),
            },
        )
    if outcome == CancelOutcome.ALREADY_DONE:
        raise HTTPException(
            status_code=409,
            detail={"code": "already_done", "message": "That job has already finished."},
        )
    return manager.record(job_id)


@router.delete("/{job_id}", status_code=204, dependencies=_WRITE_GUARDS)
def forget_job(job_id: str, manager: JobManager = Depends(get_job_manager)) -> Response:
    """Drop a finished job from the in-memory registry."""
    _job_or_404(manager, job_id)
    try:
        manager.forget(job_id)
    except ValueError:
        raise HTTPException(
            status_code=409,
            detail={
                "code": ErrorCode.CONFLICT.value,
                "message": "That job is still running.",
            },
        ) from None
    return Response(status_code=204)


# ---------------------------------------------------------------------------
# SSE plumbing
# ---------------------------------------------------------------------------


def _resume_cursor(since: int, last_event_id: str | None) -> int:
    """`Last-Event-ID` wins over `?since=`; both name the same global seq."""
    if last_event_id:
        try:
            return int(last_event_id)
        except ValueError:
            pass
    return since


def _event_payload(manager: JobManager, event: JobEvent) -> dict[str, Any]:
    """The `data:` object for one event.

    `done` is the exception: it carries the job's complete record, so a client
    that has just connected does not need a follow-up fetch to learn the
    outcome.
    """
    if event.type == "done":
        try:
            payload = manager.record(event.job_id).model_dump(mode="json")
        except KeyError:
            # The job was forgotten between the append and this read.
            payload = event.model_dump(mode="json", exclude_none=True)
        payload["seq"] = event.seq
        payload["type"] = "done"
        payload["job_id"] = event.job_id
        return payload
    # exclude_none: one JobEvent model serves six frame types, and a `log`
    # frame carrying four null progress keys would be merged into the client's
    # progress object as real data.
    return event.model_dump(mode="json", exclude_none=True)


def _frame(name: str, payload: dict[str, Any], seq: int | None = None) -> str:
    """One SSE frame. `data:` is always a single line -- json.dumps escapes
    the newlines that a log line may well contain."""
    seq = payload.get("seq") if seq is None else seq
    prefix = f"id: {seq}\n" if seq is not None else ""
    return f"{prefix}event: {name}\ndata: {json.dumps(payload, default=str)}\n\n"
