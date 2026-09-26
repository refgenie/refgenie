"""Process-wide refgenie configuration.

This package only re-exports. The environment reads, the home-directory
``mkdir`` and the ``config`` object live in :mod:`refgenie.config.env`; the
models live in :mod:`refgenie.config.settings`; the database configuration
lives in :mod:`refgenie.config.db` (imported directly, never from here, because
it depends on ``REFGENIE_HOME_PATH`` below).
"""

from refgenie.config.env import (
    REFGENIE_DB_CONFIG_PATH,
    REFGENIE_GENOME_FOLDER,
    REFGENIE_GENOME_STAGE_FOLDER,
    REFGENIE_HOME_PATH,
    REFGENIE_LOG_LEVEL,
    config,
)
from refgenie.config.settings import LogLevel, RefgenieConfig

__all__ = [
    "REFGENIE_DB_CONFIG_PATH",
    "REFGENIE_GENOME_FOLDER",
    "REFGENIE_GENOME_STAGE_FOLDER",
    "REFGENIE_HOME_PATH",
    "REFGENIE_LOG_LEVEL",
    "LogLevel",
    "RefgenieConfig",
    "config",
]
