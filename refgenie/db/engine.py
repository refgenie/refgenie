"""Creating the database engine from the database configuration file."""

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


def load_database_config(config_path: str | Path | None = None) -> DatabaseConfig:
    """
    Get the database configuration, creating a default config file if absent.

    Args:
        config_path: The path to the database configuration file.

    Returns:
        DatabaseConfig: The database configuration.
    """
    if isinstance(config_path, str):
        config_path = Path(config_path)
    cp = config_path or config.database_config_path
    if isinstance(cp, str):
        cp = Path(cp)
    if not cp.exists():
        logger.debug(
            f"Database configuration file not found at {cp}. "
            "Creating a new one with default values."
        )
        cp.parent.mkdir(parents=True, exist_ok=True)
        cp.write_text(dump_db_config(create_default_db_config()))
        logger.info(
            f"Database configuration file created at {cp}. "
            "Please edit it to set the database connection parameters."
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
