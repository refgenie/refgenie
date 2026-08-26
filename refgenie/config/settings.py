"""Models for the process-wide refgenie configuration.

Deliberately separate from ``refgenie.models``: ``refgenie.config`` is imported
by ``refgenie.logger``, which is imported by nearly everything, and
``refgenie.models`` pulls in ``refgenie.db.tables`` and therefore all of
SQLModel and SQLAlchemy. Keeping these two classes here saves ``refgenie
--help`` about 225ms.
"""

from enum import Enum
from pathlib import Path

from pydantic import BaseModel


class LogLevel(Enum):
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"


class RefgenieConfig(BaseModel):
    """
    A model for refgenie configuration, like logging level that can be set with a `REFGENIE_LOG_LEVEL` environment variable.

    The ``bridge_*`` fields configure the localhost bridge (whether the public
    SPA at an allowlisted origin may talk to a local ``refgenie dash``); they
    mirror the ``REFGENIE_BRIDGE_*`` environment variables read in
    ``refgenie/config/__init__.py``.
    """

    log_level: LogLevel
    genome_folder: Path
    genome_stage_folder: Path
    database_config_path: Path
    bridge_mode: str = "read"
    bridge_origins: str = "https://refgenie.org,https://ui.refgenie.org"
    bridge_origin_regex: str = ""
    bridge_expose_paths: bool = False
