"""The ``refgenie`` logger, and the console setup that entry points apply.

Importing this module changes no logging state. refgenie is also a library
(looper imports its pre-submit hook), and a library must leave the logging of
the program that imports it alone. Only the commands that own their process
call :func:`configure_logging`: the ``refgenie`` and ``refgenie-mcp`` commands
and the apps uvicorn loads from ``refgenie.server.main``.
"""

import logging
import sys

from refgenie.config import config

logger = logging.getLogger("refgenie")


class _StderrHandler(logging.StreamHandler):
    """Writes to whatever ``sys.stderr`` is now, so capture and redirection see it."""

    @property
    def stream(self):
        return sys.stderr

    @stream.setter
    def stream(self, value):
        pass


def configure_logging() -> None:
    """Print log records to stderr as plain lines, at ``REFGENIE_LOG_LEVEL`` (WARNING by default).

    Sets up the root logger, so call it only from an entry point. Calling it
    again replaces the handler an earlier call added.
    """
    level = getattr(logging, config.log_level.value)
    handler = _StderrHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    handler.setLevel(level)

    root = logging.getLogger()
    for old in [h for h in root.handlers if isinstance(h, _StderrHandler)]:
        root.removeHandler(old)
    root.addHandler(handler)
    root.setLevel(level)
    # The dash job log (server/jobs/logcapture.py) reads INFO records from this
    # logger; the console handler above filters them out.
    logger.setLevel(min(level, logging.INFO))
    for name in ("httpcore2", "alembic"):
        logging.getLogger(name).setLevel(logging.WARNING)
    # One INFO line per HTTP request drowns refgenie's own progress; show them only when debugging.
    logging.getLogger("httpx2").setLevel(
        logging.DEBUG if level == logging.DEBUG else logging.WARNING
    )
