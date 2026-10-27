"""Route the `refgenie` logger's output to whichever job produced it.

Every manager in the package narrates its work to one logger ("Extracting
asset tarball...", "Rolling back...", "Staging asset..."). Those lines are
exactly what a user watching a pull wants to see, and none of the managers
know anything about jobs. So instead of threading a reporter through them, one
handler on the `refgenie` logger looks at which thread emitted each record and
attributes it to the job running on that thread.

Two behaviors this handler must never violate:

* it must never raise into the logging call, which would kill the operation
  being logged (`logging` swallows handler errors only if `raiseExceptions`
  is off, which is not something to rely on); and
* it must never re-enter, because the callback it invokes is free to log.
"""

import logging
import threading

from refgenie.logger import logger as refgenie_logger

__all__ = ["JobLogHandler"]


class JobLogHandler(logging.Handler):
    """Attribute `refgenie` log records to the job running on the calling thread.

    Args:
        route: Callable taking (thread_id, level_name, message). It returns
            nothing and is expected to no-op when the thread owns no job.
        level: Minimum level to forward.
    """

    def __init__(self, route, level: int = logging.INFO):
        super().__init__(level=level)
        self._route = route
        self._reentrant = threading.local()

    def emit(self, record: logging.LogRecord) -> None:
        # Re-entrancy guard: `route` appends an event, and anything it touches
        # (or a future observer of it) may log. Without this, one log line can
        # recurse until the stack runs out.
        if getattr(self._reentrant, "active", False):
            return
        self._reentrant.active = True
        try:
            self._route(threading.get_ident(), record.levelname, record.getMessage())
        except Exception:
            # A broken sink must not take down the pull it is describing.
            self.handleError(record)
        finally:
            self._reentrant.active = False


def install(handler: JobLogHandler) -> None:
    """Attach `handler` to the package logger.

    Note for readers chasing a missing build line: pypiper calls the
    process-global `logging.disable()` around `stop_pipeline`, so a handful of
    build-teardown records never reach any handler. That is pypiper's choice,
    not a bug here; the full text is in the pipeline log file, which the tailer
    reads directly.
    """
    refgenie_logger.addHandler(handler)


def remove(handler: JobLogHandler) -> None:
    """Detach `handler` from the package logger. Safe if never installed."""
    refgenie_logger.removeHandler(handler)
