"""The jobs HTTP API, mounted at ``/v1`` in local mode.

``GET /v1/jobs/events`` is one SSE stream multiplexed across every job:

* ``event:`` is one of ``status``, ``progress``, ``log``, ``done``,
  ``truncated``, ``heartbeat``; ``data:`` is one line of JSON without null fields.
* ``id:`` is the manager-global ``seq``; resume with ``Last-Event-ID`` or ``?since=``.
* Each job gets exactly one ``done``, carrying its full ``JobRecord``.
* ``heartbeat`` is sent after 15 s of quiet; ``truncated`` carries ``dropped``
  when the ring lost unseen events.

See "Background jobs in local mode" in docs/design-notes.md.
"""

import asyncio
import itertools
import json
from collections.abc import Iterator
from typing import Any, TextIO

from fastapi import APIRouter, Depends, Header, HTTPException, Query, Request, Response
from fastapi.responses import StreamingResponse

from refgenie.models import Page
from refgenie.server.errors import ErrorCode
from refgenie.server.jobs.manager import CancelOutcome, JobManager
from refgenie.server.jobs.schemas import (
    JobEvent,
    JobKind,
    JobRecord,
    JobStatus,
    JobStatusFilter,
    utcnow,
)
from refgenie.server.local.security import require_action_header, require_action_origin
from refgenie.server.pagination import PageParams, PageResponse, get_page_params, respond

__all__ = ["get_job_manager", "router"]

router = APIRouter(prefix="/jobs", tags=["Jobs"])


#: Guards for the write routes only; the SSE stream cannot send headers.
_WRITE_GUARDS = [Depends(require_action_header), Depends(require_action_origin)]

#: How often the SSE generator polls the in-memory event ring.
POLL_SECONDS = 0.25
#: Seconds of quiet before a heartbeat frame.
HEARTBEAT_SECONDS = 15.0

#: Cap on characters returned per log line, so `get_job_log` memory is bounded.
_LOG_LINE_MAX_CHARS = 4096
_LOG_LINE_TRUNCATED_SUFFIX = " … [line truncated]"


def get_job_manager(request: Request) -> JobManager:
    """The app's JobManager, or a 503 in server mode, which has none."""
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
# Reads
# ---------------------------------------------------------------------------


@router.get("", response_model=PageResponse[JobRecord])
def list_jobs(
    kind: JobKind | None = None,
    status: JobStatusFilter | None = None,
    paging: PageParams = Depends(get_page_params),
    manager: JobManager = Depends(get_job_manager),
) -> PageResponse[JobRecord]:
    """Job records, newest first.

    `status` also accepts `active` (queued or running) and `terminal`
    (succeeded, failed or cancelled).
    """
    records = manager.list(kind=kind, status=_status_filter(status))
    page = Page(
        items=records[paging.offset : paging.offset + paging.page_size],
        offset=paging.offset,
        limit=paging.page_size,
        total=len(records),
    )
    return respond(page, paging, lambda record: record)


# Registered BEFORE /{job_id} so "events" is never captured as a job id.
@router.get("/events")
async def stream_events(
    request: Request,
    since: int = 0,
    last_event_id: str | None = Header(None, alias="Last-Event-ID"),
    manager: JobManager = Depends(get_job_manager),
) -> StreamingResponse:
    """The multiplexed SSE stream for all jobs; see the module docstring.

    Replays from the cursor, so a client connecting after a job finished still
    gets its `done`.
    """
    cursor = _resume_cursor(since, last_event_id)

    async def generate():
        nonlocal cursor
        quiet = 0.0
        while True:
            if await request.is_disconnected():
                return
            events, cursor, dropped = manager.events_since(cursor)
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
                yield _frame("heartbeat", {"seq": cursor, "ts": utcnow().isoformat()}, seq=cursor)

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
    """Polling fallback: this job's events since `since`, the next cursor, and its record."""
    record = _job_or_404(manager, job_id)
    events, cursor, dropped = manager.events_since(since)
    return {
        # Same serializer as the stream, so frames are identical.
        "events": [_event_payload(manager, event) for event in events if event.job_id == job_id],
        "next": cursor,
        "truncated": dropped > 0,
        "job": record.model_dump(mode="json"),
    }


def _bounded_lines(handle: TextIO, max_chars: int) -> Iterator[str]:
    """Yield each line without its newline, cut to `max_chars`.

    Reads in `max_chars` chunks so one giant line is never held in memory whole;
    the rest of an overlong line is dropped.
    """
    while True:
        chunk = handle.readline(max_chars)
        if chunk == "":
            return
        if chunk.endswith("\n"):
            yield chunk[:-1]
            continue
        # Hit the cap or EOF without a newline; drain to the real end of line.
        drained = False
        while True:
            more = handle.readline(max_chars)
            if more == "":
                break
            drained = True
            if more.endswith("\n"):
                break
        yield chunk + _LOG_LINE_TRUNCATED_SUFFIX if drained else chunk


def _log_window(handle: TextIO, offset: int, limit: int) -> tuple[list[str], bool]:
    """The `[offset, offset + limit)` window of lines, and whether more follow."""
    lines = _bounded_lines(handle, _LOG_LINE_MAX_CHARS)
    window = list(itertools.islice(lines, offset, offset + limit))
    has_more = next(lines, None) is not None
    return window, has_more


@router.get("/{job_id}/log", response_model=dict)
def get_job_log(
    job_id: str,
    offset: int = Query(0, ge=0),
    limit: int = Query(500, ge=1, le=5000),
    manager: JobManager = Depends(get_job_manager),
) -> dict[str, Any]:
    """Lines from the job's pipeline log file, starting at line index `offset`.

    Only the requested window is read; overlong lines are cut. `has_more` is
    true when lines exist after this page.
    """
    record = _job_or_404(manager, job_id)
    if not record.log_file:
        return {"lines": [], "next_offset": offset, "has_more": False}
    try:
        with open(record.log_file, "r", errors="replace") as handle:
            window, has_more = _log_window(handle, offset, limit)
    except OSError:
        return {"lines": [], "next_offset": offset, "has_more": False}
    return {"lines": window, "next_offset": offset + len(window), "has_more": has_more}


# ---------------------------------------------------------------------------
# Writes
# ---------------------------------------------------------------------------


@router.post(
    "/{job_id}/cancel", status_code=202, response_model=JobRecord, dependencies=_WRITE_GUARDS
)
def cancel_job(job_id: str, manager: JobManager = Depends(get_job_manager)) -> JobRecord:
    """Ask a job to stop.

    409 `not_cancellable` for a running build, `already_done` for a finished job.
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
    """The `data:` object for one event; `done` carries the full job record."""
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
    # Drop nulls so a `log` frame carries no empty progress keys.
    return event.model_dump(mode="json", exclude_none=True)


def _frame(name: str, payload: dict[str, Any], seq: int | None = None) -> str:
    """One SSE frame; json.dumps keeps `data:` on a single line."""
    seq = payload.get("seq") if seq is None else seq
    prefix = f"id: {seq}\n" if seq is not None else ""
    return f"{prefix}event: {name}\ndata: {json.dumps(payload, default=str)}\n\n"
