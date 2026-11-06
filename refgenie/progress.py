"""Optional progress reporting for long-running library operations.

Library code calls `emit(...)`; with no sink installed it is a no-op. A caller
that wants progress (the local web UI's job manager) installs a sink for one
operation with `use_sink`. Raising `ProgressAborted` from a sink cancels the
operation at its next report. Stdlib-only, so it works without the extras.
See "Progress reporting" in docs/design-notes.md.
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

    Not a `RefgenieError`, so `except RefgenieError` does not swallow it.
    """


_sink: ContextVar[ProgressSink | None] = ContextVar("refgenie_progress_sink", default=None)


def active_sink() -> ProgressSink | None:
    """The sink installed for this context, or None (lets callers skip a `rich` display)."""
    return _sink.get()


@contextmanager
def use_sink(sink: ProgressSink) -> Iterator[None]:
    """Install `sink` for the duration of the block.

    Always reset, even on exception, because pool threads are reused.
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
        ProgressAborted: If the sink raises it to cancel the operation. Other
            sink exceptions propagate too.
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
