"""The `asset-class` command group: models and handlers."""

from collections.abc import Callable

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliPositionalArg, CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.errors import fail


class AssetClassShowModel(BaseModel):
    """asset_class show: display an asset class."""

    asset_class_name: CliPositionalArg[str] = Field(description="Asset class name.")
    asset_class_version: str | None = Field(
        None,
        description="Asset class version.",
        validation_alias=AliasChoices("asset-class-version"),
    )


class AssetClassAddModel(BaseModel):
    """asset_class add: add an asset class from a source."""

    source: str = Field(description="Path/URL to the asset class to add.")
    force: bool = Field(
        False,
        description="Force the action.",
        validation_alias=AliasChoices("f", "force"),
    )


class AssetClassRemoveModel(BaseModel):
    """asset_class remove: remove an asset class."""

    asset_class_name: CliPositionalArg[str] = Field(description="Asset class name.")
    asset_class_version: str | None = Field(
        None,
        description="Asset class version.",
        validation_alias=AliasChoices("asset-class-version"),
    )


class AssetClassListModel(BaseModel):
    """asset_class list: list local asset classes."""

    pass


class AssetClassParser(BaseModel):
    """Intermediate parser for asset_class subcommands."""

    show: CliSubCommand[AssetClassShowModel] = Field(description="Show asset classes.")
    add: CliSubCommand[AssetClassAddModel] = Field(description="Add asset classes.")
    remove: CliSubCommand[AssetClassRemoveModel] = Field(description="Remove asset classes.")
    list: CliSubCommand[AssetClassListModel] = Field(description="List asset classes.")


def handle_asset_class_list(cmd, refgenie) -> None:
    rprint(refgenie.asset_class.table())


def handle_asset_class_add(cmd, refgenie) -> None:
    refgenie.asset_class.add(asset_class_source=cmd.source, exists_overwrite=cmd.force)


def handle_asset_class_remove(cmd, refgenie) -> None:
    refgenie.asset_class.remove(
        asset_class_name=cmd.asset_class_name,
        asset_class_version=cmd.asset_class_version,
    )


def handle_asset_class_show(cmd, refgenie) -> None:
    from refgenie.utils.io import cli_show_yaml

    cli_show_yaml(
        refgenie.asset_class.get(
            asset_class_name=cmd.asset_class_name,
            asset_class_version=cmd.asset_class_version,
        ).to_yaml()
    )


ASSET_CLASS_DISPATCH: dict[type, Callable] = {
    AssetClassListModel: handle_asset_class_list,
    AssetClassAddModel: handle_asset_class_add,
    AssetClassRemoveModel: handle_asset_class_remove,
    AssetClassShowModel: handle_asset_class_show,
}


def handle_asset_class_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = ASSET_CLASS_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown asset_class subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
