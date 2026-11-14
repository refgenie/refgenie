"""The `config` command group: models and handlers."""

import sys
from collections.abc import Callable

import yaml
from pydantic import BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand

from refgenie.cli.errors import fail
from refgenie.config import config
from refgenie.exceptions import MissingConfigDataError


class ConfigGetModel(BaseModel):
    """config get: display configuration."""

    pass


class ConfigSetModel(BaseModel):
    """config set: modify configuration (not yet implemented)."""

    pass


class ConfigParser(BaseModel):
    """Intermediate parser for config subcommands."""

    get: CliSubCommand[ConfigGetModel] = Field(description="Get config.")
    set: CliSubCommand[ConfigSetModel] = Field(description="Set config. [not implemented]")


def handle_config_get(cmd, refgenie) -> None:
    try:
        row = refgenie.configuration.get_latest()
    except MissingConfigDataError:
        fail("No configuration is set. Run 'refgenie init' first.")
    sys.stdout.write(
        yaml.safe_dump(
            {
                "configuration": row.model_dump(
                    mode="json", exclude={"id", "created_at", "updated_at"}
                ),
                "environment": config.model_dump(mode="json"),
            },
            sort_keys=False,
        )
    )


def handle_config_set(cmd, refgenie) -> None:
    fail("Setting configuration is not implemented yet.")


CONFIG_DISPATCH: dict[type, Callable] = {
    ConfigGetModel: handle_config_get,
    ConfigSetModel: handle_config_set,
}


def handle_config_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = CONFIG_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown config subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
