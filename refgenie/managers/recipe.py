from pathlib import Path
from collections.abc import Iterable
from typing import Any

from rich.table import Table
from sqlalchemy.exc import NoResultFound
from sqlalchemy.orm import selectinload
from sqlmodel import join, select

from refgenie.db.tables import (
    Asset,
    AssetClass,
    Recipe,
    RecipeAssetClassesInputs,
)
from refgenie.exceptions import (
    ConfigError,
    MissingAssetClassError,
    MissingRecipeError,
    RecipeExistsError,
)
from refgenie.logger import logger
from refgenie.models import InputEntity
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.utils.io import read_yaml
from refgenie.utils.tables import build_table
from refgenie.utils.versioning import parse_name_version, resolve_latest, validate_semver


class RecipeManager(ResourceManager):
    """
    A manager for recipes.
    """

    @staticmethod
    def _resolve_asset_class(session, name: str, version: str | None, label: str) -> AssetClass:
        """
        Look up the asset class a recipe names, at a pinned or the latest version.

        Args:
            session: The session to query in.
            name: The asset class name, with any version already split off.
            version: The pinned version, or None for the latest by semver.
            label: The spelling the recipe used, for the error message.

        Returns:
            AssetClass: The matching asset class.

        Raises:
            MissingAssetClassError: If no such asset class is registered.
        """
        if version is None:
            try:
                return resolve_latest(session, AssetClass, name)
            except NoResultFound:
                raise MissingAssetClassError(label) from None
        return one_or_raise(
            session,
            select(AssetClass).where(AssetClass.name == name, AssetClass.version == version),
            MissingAssetClassError(label),
            unique=True,
        )

    def exists(self, recipe_name: str, recipe_version: str | None = None) -> bool:
        """
        Check if a recipe exists.

        Args:
            recipe_name: The name of the recipe.
            recipe_version: The version of the recipe.

        Returns:
            bool: Whether the recipe exists.
        """
        where_clause = [Recipe.name == recipe_name]
        if recipe_version:
            where_clause.append(Recipe.version == recipe_version)
        with self._database_session as session:
            result = session.exec(select(Recipe).where(*where_clause))
            return bool(result.first())

    def add(self, recipe_source: Path | str, exists_overwrite: bool = False) -> Recipe:
        """
        Read a recipe from YAML file and register.

        Args:
            recipe_source: The path to the recipe file or URL.

        Returns:
            Recipe: The registered recipe.
        """

        logger.debug(f"Loading recipe from {recipe_source}")
        recipe_dict = read_yaml(recipe_source)
        validated_recipe = Recipe.model_validate(recipe_dict)
        logger.debug(f"Validated recipe: {validated_recipe}")
        recipe_name, recipe_version = recipe_dict["name"], recipe_dict["version"]

        if not validate_semver(recipe_version):
            raise ValueError(
                f"Invalid version '{recipe_version}' for recipe '{recipe_name}'. "
                f"Version must be valid semver (e.g., '0.1.0', '1.2.3')."
            )
        if not (output_asset_class_name := recipe_dict.get("output_asset_class")):
            raise ValueError("output_asset_class is required")

        if self.exists(recipe_name, recipe_version):
            if not exists_overwrite:
                raise RecipeExistsError(recipe_name, recipe_version)
            logger.info(f"Recipe '{recipe_name}' already exists. Overwriting.")
            self.remove(recipe_name, recipe_version)

        with self._database_session as session:
            # Parse output_asset_class — supports "name" or "name:version"
            oac_name, oac_version = parse_name_version(output_asset_class_name)
            asset_class = self._resolve_asset_class(
                session, oac_name, oac_version, output_asset_class_name
            )
            recipe = Recipe.model_validate(recipe_dict)
            session.add(recipe)
            recipe.output_asset_class = asset_class
            if (input_asset_spec := recipe_dict.get("input_assets")) is not None:
                recipe.input_asset_classes = []
                for input_asset_name, input_asset in input_asset_spec.items():
                    # Parse input asset class — supports "name" or "name:version"
                    iac_name, iac_version = parse_name_version(input_asset["asset_class"])
                    input_asset_class = self._resolve_asset_class(
                        session, iac_name, iac_version, input_asset["asset_class"]
                    )
                    recipe.input_asset_classes.append(
                        RecipeAssetClassesInputs(
                            asset_class=input_asset_class,
                            name=input_asset_name,
                            default=input_asset["default"],
                        )
                    )
                    session.add(input_asset_class)
            session.commit()
            logger.info(f"Registered '{recipe.name}' recipe")
            return recipe

    def list_by_output_asset_class(self, output_asset_class: str) -> Iterable[Recipe]:
        """
        Given an output asset class name, list all recipes that produce it.

        Returns recipes for all versions of the named asset class.

        Args:
            output_asset_class: The name of the output asset class.

        Returns:
            list[Recipe]: A list of all recipes that produce the output asset class.
        """
        with self._database_session as session:
            asset_class_ids = (
                session.exec(select(AssetClass.id).where(AssetClass.name == output_asset_class))
                .unique()
                .all()
            )
            if not asset_class_ids:
                return []
            return (
                session.exec(
                    select(Recipe)
                    .where(Recipe.output_asset_class_id.in_(asset_class_ids))
                    .options(
                        selectinload(Recipe.output_asset_class),
                        selectinload(Recipe.input_asset_classes),
                    )
                )
                .unique()
                .all()
            )

    def get(self, recipe_name: str, recipe_version: str | None = None) -> Recipe:
        """
        Given a recipe name and possibly a version, return the recipe.

        If version is omitted, returns the latest version by semver ordering.

        Args:
            recipe_name: The name of the recipe.
            recipe_version: The version of the recipe. If None, returns latest.
        """
        label = recipe_name if not recipe_version else f"{recipe_name} v{recipe_version}"
        with self._database_session as session:
            if recipe_version is not None:
                statement = select(Recipe).where(
                    Recipe.name == recipe_name, Recipe.version == recipe_version
                )
            else:
                try:
                    resolved = resolve_latest(session, Recipe, recipe_name)
                except NoResultFound:
                    raise MissingRecipeError(label) from None
                statement = select(Recipe).where(Recipe.id == resolved.id)
            return one_or_raise(
                session,
                statement.options(
                    selectinload(Recipe.output_asset_class).selectinload(AssetClass.seek_keys)
                ).options(selectinload(Recipe.input_asset_classes)),
                MissingRecipeError(label),
                unique=True,
            )

    def remove(self, recipe_name: str, recipe_version: str | None = None):
        """
        Remove a recipe.

        If version is omitted and only one version exists, removes it.
        If version is omitted and multiple versions exist, raises an error.

        Args:
            recipe_name: The name of the recipe.
            recipe_version: The version of the recipe.
        """
        assets_built_with_recipe = self.get_assets_built_with_recipe(
            recipe_name=recipe_name, recipe_version=recipe_version
        )
        if assets_built_with_recipe:
            raise ConfigError(
                f"Cannot remove recipe '{recipe_name}' because it is used to build the"
                f"following assets: {', '.join([asset.name for asset in assets_built_with_recipe])}"
            )
        with self._database_session as session:
            if recipe_version is not None:
                recipe = one_or_raise(
                    session,
                    select(Recipe).where(
                        Recipe.name == recipe_name, Recipe.version == recipe_version
                    ),
                    MissingRecipeError(f"{recipe_name} v{recipe_version}"),
                    unique=True,
                )
            else:
                results = (
                    session.exec(select(Recipe).where(Recipe.name == recipe_name)).unique().all()
                )
                if not results:
                    raise MissingRecipeError(recipe_name)
                if len(results) > 1:
                    versions = [r.version for r in results]
                    raise ConfigError(
                        f"Multiple versions of recipe '{recipe_name}' exist: "
                        f"{versions}. Specify a version to remove."
                    )
                recipe = results[0]
            logger.debug(f"Removing recipe: {recipe}")
            session.delete(recipe)
            session.commit()
        logger.info(f"Removed recipe '{recipe_name}'")

    def list_all(self) -> Iterable[Recipe]:
        """
        List all recipes.

        Returns:
            list[Recipe]: A list of all recipes.
        """
        with self._database_session as session:
            return (
                session.exec(
                    select(Recipe).options(
                        selectinload(Recipe.output_asset_class),
                        selectinload(Recipe.input_asset_classes),
                    )
                )
                .unique()
                .all()
            )

    def table(self, recipe_names: list[str] | None = None) -> Table:
        """
        Create a table of all recipes.
        """
        recipes = self.list_all()
        if recipe_names is not None:
            recipes = [recipe for recipe in recipes if recipe.name in recipe_names]
        return build_table(
            "Recipes",
            [
                "Name",
                "Version",
                "Output asset class",
                "Input asset classes",
                "Input files",
                "Input params",
                "Docker image",
            ],
            [
                (
                    recipe.name,
                    recipe.version,
                    recipe.output_asset_class.name if recipe.output_asset_class else "N/A",
                    self.input_enities_to_text(recipe.input_assets),
                    self.input_enities_to_text(recipe.input_files),
                    self.input_enities_to_text(recipe.input_params),
                    recipe.docker_image,
                )
                for recipe in recipes
            ],
            end_section=True,
        )

    def get_required_input_asset_classes(
        self,
        recipe_name: str,
        recipe_version: str | None = None,
    ) -> Iterable[tuple[AssetClass, str]]:
        """
        Given a recipe name list all input asset classes and possibly their defaults

        Args:
            recipe_name: The name of the recipe.
            recipe_version: The version of the recipe.

        Returns:
            list of (AssetClass, default) tuples for the recipe's inputs.
        """
        recipe = self.get(recipe_name=recipe_name, recipe_version=recipe_version)
        with self._database_session as session:
            input_asset_classes = (
                session.exec(
                    select(AssetClass, RecipeAssetClassesInputs.default)
                    .select_from(join(RecipeAssetClassesInputs, AssetClass))
                    .where(RecipeAssetClassesInputs.recipe_id == recipe.id)
                )
                .unique()
                .all()
            )
            return input_asset_classes

    def get_assets_built_with_recipe(
        self, recipe_name: str, recipe_version: str | None = None
    ) -> Iterable[Asset]:
        """
        Get all assets that are built with the recipe.

        Args:
            recipe_name: The name of the recipe.
            recipe_version: The version of the recipe.

        Returns:
            list[Asset]: The list of assets that are built with the recipe.
        """
        recipe = self.get(recipe_name=recipe_name, recipe_version=recipe_version)
        with self._database_session as session:
            return session.exec(select(Asset).where(Asset.recipe_id == recipe.id)).unique().all()

    @staticmethod
    def input_enities_to_text(
        input_entities: dict[str, dict[str, Any | None]] = None,
    ) -> str:
        def _make_bullet_list(items: list[str]) -> str:
            details_text = []
            for item in items:
                details_text.append(f"• {item}")
            return "\n".join(details_text)

        if not input_entities:
            return "[dim]None[/dim]"
        strings = []
        for input_entity_id, input_entity_data in input_entities.items():
            strings.append(f"{input_entity_id} [dim]{InputEntity(**input_entity_data)}[/dim]")
        return _make_bullet_list(strings)
