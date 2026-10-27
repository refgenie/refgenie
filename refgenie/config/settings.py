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
    """

    log_level: LogLevel
    genome_folder: Path
    genome_stage_folder: Path
    database_config_path: Path
    #: `REFGENIE_DISABLE_PLUGINS`: "1"/"true"/"all" turns off every plugin; a
    #: comma-separated list of plugin names turns off those.
    disable_plugins: str = ""
    #: `REFGENIE_SERVER_PLUGINS`: run plugins in server mode too.
    server_plugins: bool = False
