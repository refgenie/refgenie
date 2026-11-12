"""Logging setup, run once at import: plain lines on stderr, WARNING by default."""

import logging
import sys

from refgenie.config import config


class _StderrHandler(logging.StreamHandler):
    """Writes to whatever ``sys.stderr`` is now, so capture and redirection see it."""

    @property
    def stream(self):
        return sys.stderr

    @stream.setter
    def stream(self, value):
        pass


_level = getattr(logging, config.log_level.value)
_handler = _StderrHandler()
_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
_handler.setLevel(_level)

logging.getLogger().handlers.clear()
logging.getLogger().addHandler(_handler)
logging.getLogger().setLevel(_level)
# The dash job log (server/jobs/logcapture.py) reads INFO records from this
# logger; the console handler above filters them out.
logging.getLogger("refgenie").setLevel(min(_level, logging.INFO))
for _name in ("httpx", "httpcore", "alembic"):
    logging.getLogger(_name).setLevel(logging.WARNING)

logger = logging.getLogger("refgenie")
