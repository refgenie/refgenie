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


class ConfigParser(BaseModel):
    """Show refgenie's configuration (`refgenie config get`, or plain `refgenie config`)."""

    get: CliSubCommand[ConfigGetModel] = Field(description="Get config.")


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


CONFIG_DISPATCH: dict[type, Callable] = {
    ConfigGetModel: handle_config_get,
}


def handle_config_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=False)
    if leaf is None:  # bare `refgenie config` prints it, like `config get`
        leaf = ConfigGetModel()
    handler = CONFIG_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown config subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
