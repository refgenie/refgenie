from enum import StrEnum
from pathlib import Path
from typing import Annotated, Literal

from pydantic import AfterValidator, BaseModel, TypeAdapter, computed_field
import yacman

from yaml import dump as ydump

from refgenie.config import REFGENIE_HOME_PATH


def expand_and_check(value: str) -> str:
    """
    Expand user directory and resolve to absolute path.
    """
    expanded_path = Path(value).expanduser().resolve()
    if not expanded_path.parent.exists() or not expanded_path.parent.is_dir():
        raise ValueError(f"Path is invalid: {expanded_path}. Parent directory does not exist.")
    return expanded_path.as_posix()


class DatabaseType(StrEnum):
    """
    Definition of supported database types.
    """

    SQLITE = "sqlite"
    POSTGRESQL = "postgresql"


class SQLiteConfig(BaseModel):
    type: Literal[DatabaseType.SQLITE] = DatabaseType.SQLITE
    path: Annotated[str, AfterValidator(expand_and_check)] = (
        REFGENIE_HOME_PATH / "refgenie"
    ).as_posix()

    @computed_field  # type: ignore[prop-decorator]
    @property
    def url(self) -> str:
        """
        The URL of the SQLite database.
        """
        return f"sqlite:///{self.path}"


class PostgresConfig(BaseModel):
    type: Literal[DatabaseType.POSTGRESQL] = DatabaseType.POSTGRESQL
    name: str
    user: str
    password: str
    host: str
    port: int

    @computed_field  # type: ignore[prop-decorator]
    @property
    def url(self) -> str:
        """
        The URL of the PostgreSQL database.
        """
        return (
            f"postgresql+psycopg://{self.user}:{self.password}@{self.host}:{self.port}/{self.name}"
        )


DatabaseConfig = Annotated[SQLiteConfig | PostgresConfig, "type"]


def load_db_config(path: Path) -> DatabaseConfig:
    """
    Load the database configuration from a YAML file using discriminated union.

    Args:
        path: The path to the YAML file.

    Returns:
        DatabaseConfig: The database configuration.
    """

    data = yacman.YAMLConfigManager.from_yaml_file(path.as_posix()).exp

    # Every SQLiteConfig field has a default, so any stray YAML (e.g. a legacy
    # genome config) would otherwise validate as an empty default DB. Require
    # the file to actually declare a database, so a wrong REFGENIE_DB_CONFIG_PATH
    # fails loudly instead of silently opening a fresh empty database.
    if not isinstance(data, dict) or "type" not in data:
        raise ValueError(
            f"Not a valid refgenie database config: {path}. "
            f"Expected a 'type' key ('sqlite' or 'postgresql')."
        )
    if data["type"] == DatabaseType.SQLITE and "path" not in data:
        raise ValueError(
            f"Invalid sqlite database config: {path}. Expected a 'path' key."
        )

    adapter: TypeAdapter = TypeAdapter(DatabaseConfig)
    return adapter.validate_python(data)


def dump_db_config(config: DatabaseConfig) -> str:
    """
    Serialize the database configuration to a YAML string.

    Args:
        config: The database configuration to save.

    Returns:
        str: The YAML representation of the database configuration.
    """
    return ydump(
        config.model_dump(mode="json", exclude_none=True, exclude={"url"}),
        default_flow_style=False,
    )


def create_default_db_config() -> DatabaseConfig:
    """
    Create a default database configuration.

    Returns:
        DatabaseConfig: The default database configuration.
    """
    return TypeAdapter(DatabaseConfig).validate_python(SQLiteConfig())
