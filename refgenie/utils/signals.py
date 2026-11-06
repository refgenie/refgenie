"""Ctrl-C handling for long-running work."""

import signal
import sys
import threading
from collections.abc import Iterator
from contextlib import contextmanager

from refgenie.logger import logger


@contextmanager
def sigint_guard(message: str) -> Iterator[None]:
    """
    On Ctrl-C inside the block, log ``message`` and exit with status 130.

    ``signal.signal`` works only on the main thread, so worker threads (web UI
    jobs) run the block unguarded. The previous handler is restored on exit.
    """
    if threading.current_thread() is not threading.main_thread():
        yield
        return

    def handle(sig, frame):
        logger.warning(f"\n{message}")
        sys.exit(130)

    previous = signal.signal(signal.SIGINT, handle)
    try:
        yield
    finally:
        signal.signal(signal.SIGINT, previous)
