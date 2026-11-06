"""Models for the process-wide refgenie configuration.

Kept out of ``refgenie.models``, which imports SQLModel, because nearly every
module imports this one through ``refgenie.logger``.
"""

from enum import Enum
from pathlib import Path
from typing import Literal

from pydantic import BaseModel


class LogLevel(Enum):
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"


class RefgenieConfig(BaseModel):
    """Process-wide settings, read from ``REFGENIE_*`` environment variables."""

    log_level: LogLevel
    genome_folder: Path
    genome_stage_folder: Path
    database_config_path: Path
    #: `REFGENIE_DISABLE_PLUGINS`: "1"/"true"/"all" turns off every plugin; a
    #: comma-separated list of plugin names turns off those.
    disable_plugins: str = ""
    #: `REFGENIE_SERVER_PLUGINS`: run plugins in server mode too.
    server_plugins: bool = False
    #: `REFGENIE_SUBSCRIBE_DEFAULT_SERVER`: subscribe a newly created
    #: configuration to the default server (https://api.refgenie.org).
    subscribe_default_server: bool = True
    #: `REFGENIE_CONTAINER_ENGINE`: what runs a recipe's container image.
    #: `auto` takes Docker when its daemon answers, else Apptainer.
    container_engine: Literal["auto", "docker", "apptainer"] = "auto"
    #: `REFGENIE_APPTAINER_CACHEDIR`: passed to Apptainer as `APPTAINER_CACHEDIR`,
    #: so image pulls land on scratch rather than in `$HOME`.
    apptainer_cachedir: Path | None = None
