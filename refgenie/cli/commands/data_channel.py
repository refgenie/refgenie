"""The `data-channel` command group: models and handlers."""

import sys
from collections.abc import Callable
from typing import Literal

from pydantic import AliasChoices, BaseModel, ConfigDict, Field, model_validator
from pydantic_settings import CliPositionalArg, CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.errors import EXIT_NOT_FOUND, fail
from refgenie.logger import logger


class DataChannelAddModel(BaseModel):
    """data_channel add: add a data channel."""

    name: CliPositionalArg[str] = Field(description="Name of the data channel.")
    type: CliPositionalArg[Literal["ftp", "http", "https", "local"]] = Field(
        description="Type of the data channel.",
    )
    index_address: CliPositionalArg[str] = Field(
        description="Address of the data channel index YAML file.",
    )
    description: str | None = Field(
        None,
        description="Description of the data channel.",
        validation_alias=AliasChoices("d", "description"),
    )
    username: str | None = Field(None, description="Username for authentication.")
    password: str | None = Field(None, description="Password for authentication.")
    token: str | None = Field(None, description="Authentication token.")


class DataChannelRemoveModel(BaseModel):
    """data_channel remove: remove a data channel."""

    name: CliPositionalArg[str] = Field(description="Name of the data channel to remove.")


class DataChannelListModel(BaseModel):
    """data_channel list: list all data channels."""

    pass


class DataChannelShowModel(BaseModel):
    """data_channel show: show data channel details."""

    name: CliPositionalArg[str] = Field(description="Name of the data channel to show details for.")


class DataChannelValidateModel(BaseModel):
    """data_channel validate: validate a data channel."""

    name: CliPositionalArg[str] = Field(description="Name of the data channel to validate.")


class DataChannelSyncNestedModel(BaseModel):
    """data_channel sync: sync from a data channel (mutually exclusive exists flags)."""

    # See AliasGetNestedModel: required so exists_ok/exists_overwrite can be set
    # by field name rather than being silently dropped.
    model_config = ConfigDict(populate_by_name=True)

    name: CliPositionalArg[str] = Field(description="Name of the data channel to sync from.")
    exists_ok: bool = Field(
        False,
        description="Skip existing assets/recipes without error.",
        validation_alias=AliasChoices("exists-ok"),
    )
    exists_overwrite: bool = Field(
        False,
        description="Delete conflicting items before adding.",
        validation_alias=AliasChoices("exists-overwrite"),
    )

    @model_validator(mode="after")
    def check_exists_exclusion(self):
        if self.exists_ok and self.exists_overwrite:
            raise ValueError("--exists-ok and --exists-overwrite are mutually exclusive")
        return self


class DataChannelParser(BaseModel):
    """Intermediate parser for data_channel subcommands."""

    add: CliSubCommand[DataChannelAddModel] = Field(description="Add a data channel.")
    remove: CliSubCommand[DataChannelRemoveModel] = Field(description="Remove a data channel.")
    list: CliSubCommand[DataChannelListModel] = Field(description="List all data channels.")
    show: CliSubCommand[DataChannelShowModel] = Field(description="Show data channel details.")
    validate_channel: CliSubCommand[DataChannelValidateModel] = Field(
        description="Validate a data channel.",
        alias="validate",
    )
    sync: CliSubCommand[DataChannelSyncNestedModel] = Field(description="Sync from a data channel.")


def handle_data_channel_list(cmd, refgenie) -> None:
    rprint(refgenie.sources.channels_table())


def handle_data_channel_show(cmd, refgenie) -> None:
    if channel := refgenie.sources.get_channel(cmd.name):
        rprint(refgenie.sources.channels_table(channel_names=[channel.name]))
    else:
        fail(f"Data channel '{cmd.name}' not found", EXIT_NOT_FOUND)


def handle_data_channel_add(cmd, refgenie) -> None:
    from refgenie.db.tables import DataChannelType

    credentials = {}
    if cmd.username:
        credentials["username"] = cmd.username
    if cmd.password:
        credentials["password"] = cmd.password
    if cmd.token:
        credentials["token"] = cmd.token
    try:
        refgenie.sources.add_channel(
            name=cmd.name,
            type=DataChannelType(cmd.type),
            index_address=cmd.index_address,
            description=cmd.description,
            credentials=credentials if credentials else None,
        )
    except (ValueError, RuntimeError, OSError) as e:
        fail(f"Failed to add data channel: {e}")
    logger.info(f"Added data channel: {cmd.name}")


def handle_data_channel_remove(cmd, refgenie) -> None:
    if not refgenie.sources.remove_channel(cmd.name):
        fail(f"Data channel '{cmd.name}' not found", EXIT_NOT_FOUND)
    logger.info(f"Removed data channel: {cmd.name}")


def handle_data_channel_validate(cmd, refgenie) -> None:
    if not refgenie.sources.test_channel(cmd.name):
        logger.info(
            "Please check whether the resource is accessible and/or validity of the index file"
        )
        fail(f"Data channel '{cmd.name}' is invalid")
    logger.info(f"Data channel '{cmd.name}' is valid")


def handle_data_channel_sync(cmd, refgenie) -> None:
    from refgenie.exceptions import AssetClassExistsError, RecipeExistsError, RefgenieError

    if not refgenie.sources.test_channel(cmd.name):
        fail(f"Data channel '{cmd.name}' is not accessible")

    success = True
    for asset_class_url in refgenie.sources.iter_asset_classes(cmd.name):
        try:
            refgenie.asset_class.add(asset_class_url, exists_overwrite=cmd.exists_overwrite)
            logger.debug(f"Added asset class from {asset_class_url}")
        except AssetClassExistsError as e:
            if cmd.exists_ok:
                logger.warning(f"Skipping adding existing asset class from {asset_class_url}")
                continue
            else:
                logger.error(f"Failed to add asset class from {asset_class_url}: {e}")
                success = False
        except (RefgenieError, OSError, ValueError) as e:
            logger.error(f"Error adding asset class from {asset_class_url}: {e}")
            success = False

    if success:
        for recipe_url in refgenie.sources.iter_recipes(cmd.name):
            try:
                refgenie.recipe.add(recipe_url, exists_overwrite=cmd.exists_overwrite)
                logger.debug(f"Added recipe from {recipe_url}")
            except RecipeExistsError as e:
                if cmd.exists_ok:
                    logger.warning(f"Skipping adding existing recipe from {recipe_url}")
                    continue
                else:
                    logger.error(f"Failed to add recipe from {recipe_url}: {e}")
                    success = False
            except (RefgenieError, OSError, ValueError) as e:
                logger.error(f"Error adding recipe from {recipe_url}: {e}")
                success = False

    if not success:
        fail(f"Failed to sync from channel '{cmd.name}'")
    logger.info(f"Successfully synced from channel '{cmd.name}'")


DATA_CHANNEL_DISPATCH: dict[type, Callable] = {
    DataChannelListModel: handle_data_channel_list,
    DataChannelShowModel: handle_data_channel_show,
    DataChannelAddModel: handle_data_channel_add,
    DataChannelRemoveModel: handle_data_channel_remove,
    DataChannelValidateModel: handle_data_channel_validate,
    DataChannelSyncNestedModel: handle_data_channel_sync,
}


def handle_data_channel_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = DATA_CHANNEL_DISPATCH.get(type(leaf))
    if handler is None:
        logger.error(f"Unknown data_channel subcommand: {type(leaf).__name__}")
        sys.exit(1)
    handler(leaf, refgenie)
