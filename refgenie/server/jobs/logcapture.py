"""Route the `refgenie` logger's output to the job running on the emitting thread.

The handler never raises into the logging call and never re-enters.
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
        # `route` may itself log; guard against recursion.
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

    pypiper calls `logging.disable()` around `stop_pipeline`, so a few teardown
    records never arrive; they are still in the pipeline log file.
    """
    refgenie_logger.addHandler(handler)


def remove(handler: JobLogHandler) -> None:
    """Detach `handler` from the package logger. Safe if never installed."""
    refgenie_logger.removeHandler(handler)
