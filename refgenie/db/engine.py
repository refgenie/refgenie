"""Creating the database engine from the database configuration file."""

import os
import re
from pathlib import Path

from sqlalchemy.engine import Engine
from sqlmodel import create_engine

from refgenie.config import config
from refgenie.config.db import (
    DatabaseConfig,
    create_default_db_config,
    dump_db_config,
    load_db_config,
)
from refgenie.logger import logger

__all__ = ["Engine", "create_database_engine", "load_database_config"]

#: A ``$VAR`` or ``${VAR}`` that ``os.path.expandvars`` left as written, because
#: the variable is unset. The same pattern ``expandvars`` itself looks for.
_UNEXPANDED_VAR = re.compile(r"\$(\w+|\{[^}]*\})")


def _given_config_path(config_path: str | Path | None) -> Path | None:
    """The configuration file a caller named, or None if it named none.

    Environment variables in the path are expanded. Empty text, and text that
    still holds an unset variable, count as no path: looper passes
    ``${REFGENIE_DB_CONFIG_PATH}`` through as written when that variable is
    unset, and a file must never be created at that literal path.
    """
    if config_path is None:
        return None
    text = os.path.expandvars(str(config_path)).strip()
    if not text or _UNEXPANDED_VAR.search(text):
        return None
    return Path(text)


def load_database_config(config_path: str | Path | None = None) -> DatabaseConfig:
    """
    Get the database configuration, creating a default config file if absent.

    Args:
        config_path: The path to the database configuration file. None, empty
            text, or a path naming an unset environment variable means the
            configured path: ``$REFGENIE_DB_CONFIG_PATH``, else
            ``$REFGENIE_HOME_PATH/refgenie_db_config.yaml``.

    Returns:
        DatabaseConfig: The database configuration.
    """
    cp = _given_config_path(config_path) or config.database_config_path
    if not cp.exists():
        logger.debug(
            f"Database configuration file not found at {cp}. "
            "Creating a new one with default values."
        )
        cp.parent.mkdir(parents=True, exist_ok=True)
        cp.write_text(dump_db_config(create_default_db_config()))
        logger.info(
            f"Database configuration file created at {cp}. "
            "It uses SQLite by default; edit it only to switch to PostgreSQL."
        )
    return load_db_config(cp)


def create_database_engine(
    database_config_path: str | Path | None = None, echo: bool = False
) -> Engine:
    """
    Get the default database engine.

    Args:
        database_config_path: The database configuration file; defaults to the
            configured path.
        echo: Whether to echo SQL statements.

    Returns:
        Engine: The database engine.
    """
    database_config = load_database_config(config_path=database_config_path)
    # pre_ping + recycle avoid hangs on pooled connections a remote Postgres
    # silently dropped while idle. No-ops for SQLite.
    return create_engine(url=database_config.url, echo=echo, pool_pre_ping=True, pool_recycle=300)
