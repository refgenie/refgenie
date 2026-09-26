"""The `config` command group: models and handlers."""

from collections.abc import Callable

from pydantic import BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.errors import fail
from refgenie.logger import logger


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
    if (table := refgenie.configuration.table()) is not None:
        rprint(table)
    else:
        logger.warning("No configuration is set.")
        logger.info("Use 'refgenie init' to initialize the configuration.")


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
