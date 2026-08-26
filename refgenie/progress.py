"""Optional progress reporting for long-running library operations.

Library code calls `emit(...)`. With no sink installed it is a no-op, so the
CLI and every plain library consumer behave exactly as before. A caller that
wants progress -- the local web UI's job manager -- installs a sink for the
duration of one operation with `use_sink`.

Why a `ContextVar` rather than a `progress=` parameter: the sink has to reach
`RefgenieserverClient.download_with_progress` from `Refgenie.pull`, through
`AssetManager.pull` -> `AssetPuller.pull` -> `_pull_archive_mode`. Threading a
parameter through all of them is four signature changes on load-bearing public
methods for one optional feature. A ContextVar set inside a worker thread is
naturally per-thread (each thread starts with a fresh top-level context), so
two concurrent jobs get independent sinks with no locking.

One documented limitation: threads spawned *by* library internals do not
inherit the installed sink. pypiper's subprocess tee threads are the case that
matters, which is why build output is collected by tailing pypiper's log file
rather than through this mechanism.

Sinks may raise. `emit` deliberately does not swallow sink exceptions: raising
`ProgressAborted` from a sink is how a caller cancels a running operation
cooperatively at the next reporting point.

This module is stdlib-only on purpose. It is imported by library code that must
work without the ``dash`` or ``server`` extras installed.
"""

from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from typing import Any, Literal

__all__ = [
    "EventType",
    "ProgressAborted",
    "ProgressEvent",
    "ProgressSink",
    "active_sink",
    "emit",
    "use_sink",
]

#: ``stage`` marks a coarse phase change ("Downloading rCRSd/fasta"),
#: ``progress`` carries a current/total pair, ``log`` a free-text line.
EventType = Literal["stage", "progress", "log"]


@dataclass(frozen=True, slots=True)
class ProgressEvent:
    """One report from a long-running operation."""

    type: EventType
    message: str | None = None
    current: int | None = None
    total: int | None = None
    unit: str | None = None
    extra: Mapping[str, Any] = field(default_factory=dict)


ProgressSink = Callable[[ProgressEvent], None]


class ProgressAborted(Exception):
    """Raised by a sink to abort the operation reporting to it.

    Not a `RefgenieError`: it is a control-flow signal from the caller, not a
    failure of the operation, and `except RefgenieError` must not swallow it.
    """


_sink: ContextVar[ProgressSink | None] = ContextVar("refgenie_progress_sink", default=None)


def active_sink() -> ProgressSink | None:
    """The sink installed for this context, or None.

    Call sites use this to choose a reporting strategy -- notably, to skip
    building a `rich` live display when something else is consuming progress.
    """
    return _sink.get()


@contextmanager
def use_sink(sink: ProgressSink) -> Iterator[None]:
    """Install `sink` for the duration of the block.

    The token is always reset, including on exception. That matters: worker
    threads in a pool are reused, and a leaked sink would send one job's
    progress to a previous job's record.
    """
    token = _sink.set(sink)
    try:
        yield
    finally:
        _sink.reset(token)


def emit(
    type: EventType,
    message: str | None = None,
    *,
    current: int | None = None,
    total: int | None = None,
    unit: str | None = None,
    **extra: Any,
) -> None:
    """Report progress to the installed sink, if there is one.

    Args:
        type: "stage", "progress" or "log".
        message: Human-readable description.
        current: Units completed so far.
        total: Units expected in total, or None when unknown (an
            indeterminate progress bar).
        unit: What `current`/`total` count ("bytes", "steps", ...).
        **extra: Arbitrary additional fields, passed through to the sink.

    Raises:
        ProgressAborted: If the sink raises it to cancel the operation. Sink
            exceptions are deliberately not swallowed.
    """
    sink = _sink.get()
    if sink is None:
        return
    sink(
        ProgressEvent(
            type=type,
            message=message,
            current=current,
            total=total,
            unit=unit,
            extra=extra,
        )
    )
