"""The `recipe` command group: models and handlers."""

import sys
from collections.abc import Callable

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliPositionalArg, CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.errors import fail


class RecipeShowModel(BaseModel):
    """recipe show: display a recipe."""

    recipe_name: CliPositionalArg[str] = Field(description="Recipe name to perform the action on.")
    recipe_version: str | None = Field(
        None,
        description="Recipe version to perform the action on.",
        validation_alias=AliasChoices("recipe-version"),
    )


class RecipeAddModel(BaseModel):
    """recipe add: add a recipe from a source."""

    source: str = Field(description="Path/URL to the recipe to add.")
    force: bool = Field(
        False,
        description="Whether to force the action.",
        validation_alias=AliasChoices("f", "force"),
    )


class RecipeRemoveModel(BaseModel):
    """recipe remove: remove a recipe."""

    recipe_name: CliPositionalArg[str] = Field(description="Recipe name.")
    recipe_version: str | None = Field(
        None,
        description="Recipe version.",
        validation_alias=AliasChoices("recipe-version"),
    )


class RecipeListModel(BaseModel):
    """recipe list: list local recipes."""

    pass


class RecipeRequirementsModel(BaseModel):
    """recipe requirements: show recipe requirements."""

    recipe_name: CliPositionalArg[str] = Field(description="Recipe name.")
    recipe_version: str | None = Field(
        None,
        description="Recipe version.",
        validation_alias=AliasChoices("recipe-version"),
    )


class RecipeParser(BaseModel):
    """Intermediate parser for recipe subcommands."""

    show: CliSubCommand[RecipeShowModel] = Field(description="Show recipes.")
    add: CliSubCommand[RecipeAddModel] = Field(description="Add recipes.")
    remove: CliSubCommand[RecipeRemoveModel] = Field(description="Remove recipes.")
    list: CliSubCommand[RecipeListModel] = Field(description="List recipes.")
    requirements: CliSubCommand[RecipeRequirementsModel] = Field(
        description="Show recipe requirements."
    )


def handle_recipe_list(cmd, refgenie) -> None:
    rprint(refgenie.recipe.table())


def handle_recipe_add(cmd, refgenie) -> None:
    refgenie.recipe.add(recipe_source=cmd.source, exists_overwrite=cmd.force)


def handle_recipe_show(cmd, refgenie) -> None:
    from refgenie.utils.io import cli_show_yaml

    cli_show_yaml(
        refgenie.recipe.get(
            recipe_name=cmd.recipe_name, recipe_version=cmd.recipe_version
        ).to_yaml()
    )


def handle_recipe_remove(cmd, refgenie) -> None:
    refgenie.recipe.remove(recipe_name=cmd.recipe_name, recipe_version=cmd.recipe_version)


def handle_recipe_requirements(cmd, refgenie) -> None:
    recipe = refgenie.recipe.get(
        recipe_name=cmd.recipe_name, recipe_version=cmd.recipe_version
    )
    rprint(refgenie.recipe.table(recipe_names=[recipe.name]))
    sys.exit(0)


RECIPE_DISPATCH: dict[type, Callable] = {
    RecipeListModel: handle_recipe_list,
    RecipeAddModel: handle_recipe_add,
    RecipeShowModel: handle_recipe_show,
    RecipeRemoveModel: handle_recipe_remove,
    RecipeRequirementsModel: handle_recipe_requirements,
}


def handle_recipe_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = RECIPE_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown recipe subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
