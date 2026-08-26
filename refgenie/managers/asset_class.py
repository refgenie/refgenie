from pathlib import Path
from collections.abc import Iterable

from rich.table import Table
from sqlalchemy.exc import NoResultFound
from sqlalchemy.orm import selectinload
from sqlmodel import select

from refgenie.db.tables import (
    AssetClass,
    AssetClassSeekKey,
    Recipe,
    RecipeAssetClassesInputs,
)
from refgenie.exceptions import (
    AssetClassExistsError,
    ConfigError,
    MissingAssetClassError,
)
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.utils.io import read_yaml
from refgenie.utils.tables import build_table
from refgenie.utils.versioning import resolve_latest, validate_semver


class AssetClassManager(ResourceManager):
    """
    A manager for asset classes.
    """

    def add(
        self, asset_class_source: Path | str, exists_overwrite: bool = False
    ) -> AssetClass:
        """
        Read a asset class from YAML file and register.

        Args:
            asset_class_source: The path to the asset class file or URL.

        Returns:
            AssetClass: The registered asset class.
        """
        logger.debug(f"Loading asset_class from {asset_class_source}")
        asset_class_dict = read_yaml(asset_class_source)

        # Detect if the user accidentally passed a recipe YAML instead of an
        # asset class YAML.  Recipes contain fields like 'command_templates'
        # and 'output_asset_class' that asset classes never have.
        recipe_only_fields = {"command_templates", "output_asset_class"}
        found = recipe_only_fields & set(asset_class_dict)
        if found:
            raise ValueError(
                f"The provided YAML appears to be a recipe, not an asset class "
                f"(contains recipe-specific fields: {', '.join(sorted(found))}). "
                f"Use 'refgenie recipe add' to register recipes, or provide an "
                f"asset class YAML file to 'asset-class add'."
            )

        asset_class = AssetClass.model_validate(asset_class_dict)

        if not validate_semver(asset_class.version):
            raise ValueError(
                f"Invalid version '{asset_class.version}' for asset class '{asset_class.name}'. "
                f"Version must be valid semver (e.g., '0.1.0', '1.2.3')."
            )

        if self.exists(asset_class.name, asset_class.version):
            if not exists_overwrite:
                raise AssetClassExistsError(asset_class.name, asset_class.version)
            logger.info(f"Asset class '{asset_class.name}' already exists. Overwriting.")
            self.remove(asset_class.name, asset_class.version)

        seek_keys = []
        for seek_key_name, seek_key in asset_class_dict.get("seek_keys", {}).items():
            seek_keys.append(
                AssetClassSeekKey(
                    name=seek_key_name,
                    value=seek_key.get("value"),
                    description=seek_key.get("description"),
                    type=seek_key["type"],
                )
            )
        asset_class.seek_keys = seek_keys
        with self._database_session as session:
            session.add(asset_class)
            session.commit()
            session.refresh(asset_class)
            logger.info(f"Registered '{asset_class.name}' asset class")
            return asset_class

    def get(self, asset_class_name: str, asset_class_version: str | None = None) -> AssetClass:
        """
        Get an asset class by its name and optionally version.

        If version is omitted, returns the latest version by semver ordering.

        Args:
            asset_class_name: The name of the asset class.
            asset_class_version: The version of the asset class. If None, returns latest.

        Returns:
            AssetClass: The asset class.
        """
        label = (
            f"{asset_class_name} v{asset_class_version}"
            if asset_class_version
            else asset_class_name
        )
        with self._database_session as session:
            if asset_class_version is not None:
                statement = select(AssetClass).where(
                    AssetClass.name == asset_class_name,
                    AssetClass.version == asset_class_version,
                )
            else:
                try:
                    ac = resolve_latest(session, AssetClass, asset_class_name)
                except NoResultFound:
                    raise MissingAssetClassError(label) from None
                statement = select(AssetClass).where(AssetClass.id == ac.id)
            # Eagerly load seek_keys
            return one_or_raise(
                session,
                statement.options(selectinload(AssetClass.seek_keys)),
                MissingAssetClassError(label),
                unique=True,
            )

    def exists(self, asset_class_name: str, asset_class_version: str) -> bool:
        """
        Check if an asset class exists.

        Args:
            asset_class_name: The name of the asset class.
            asset_class_version: The version of the asset class.

        Returns:
            bool: Whether the asset class exists.
        """
        with self._database_session as session:
            return bool(
                session.exec(
                    select(AssetClass).where(
                        AssetClass.name == asset_class_name,
                        AssetClass.version == asset_class_version,
                    )
                ).first()
            )

    def remove(self, asset_class_name: str, asset_class_version: str | None = None):
        """
        Remove an asset class.

        Note that asset classes can be removed only if there are currently no assets of that
        class managed by refgenie, and no recipes that have inputs or outputs of that class.

        If version is omitted and only one version exists, removes it.
        If version is omitted and multiple versions exist, raises an error.

        Args:
            asset_class_name: The name of the asset class.
            asset_class_version: The version of the asset class.
        """
        logger.debug(f"Removing asset class '{asset_class_name}'")
        with self._database_session as session:
            if asset_class_version is not None:
                query = select(AssetClass).where(
                    AssetClass.name == asset_class_name,
                    AssetClass.version == asset_class_version,
                )
                asset_class = one_or_raise(
                    session,
                    query,
                    MissingAssetClassError(f"{asset_class_name} v{asset_class_version}"),
                    unique=True,
                )
            else:
                results = (
                    session.exec(select(AssetClass).where(AssetClass.name == asset_class_name))
                    .unique()
                    .all()
                )
                if not results:
                    raise MissingAssetClassError(asset_class_name)
                if len(results) > 1:
                    versions = [r.version for r in results]
                    raise ConfigError(
                        f"Multiple versions of asset class '{asset_class_name}' exist: "
                        f"{versions}. Specify a version to remove."
                    )
                asset_class = results[0]
            if asset_class.asset_groups:
                raise ConfigError(
                    f"Asset class '{asset_class_name} v{asset_class_version}' is in use by the "
                    f"following asset groups and cannot be removed: {asset_class.asset_groups}"
                )
            recipe = (
                session.exec(
                    select(Recipe).where(Recipe.output_asset_class_id == asset_class.id)
                ).first()
                or session.exec(
                    select(Recipe)
                    .join(RecipeAssetClassesInputs)
                    .where(RecipeAssetClassesInputs.asset_class_id == asset_class.id)
                ).first()
            )
            logger.debug(f"Recipe: {recipe}")
            if recipe:
                raise ConfigError(
                    f"Asset class '{asset_class_name} v{asset_class_version}' is in use by the "
                    f"following recipe and cannot be removed: '{recipe.name} v{recipe.version}'"
                )
            session.delete(asset_class)
            session.commit()
        logger.info(f"Removed asset class '{asset_class_name}'")

    def list_all(self) -> Iterable[AssetClass]:
        """
        List all asset classes.

        Returns:
            list[AssetClass]: A list of all asset classes.
        """
        with self._database_session as session:
            return (
                session.exec(
                    select(AssetClass).options(
                        selectinload(AssetClass.seek_keys),
                    )
                )
                .unique()
                .all()
            )

    def list_root(self) -> Iterable[AssetClass]:
        """
        List all root asset classes (asset classes that have no dependencies).

        Returns:
            list[AssetClass]: A list of all root asset classes.
        """
        from sqlalchemy.sql.expression import null

        with self._database_session as session:
            # asset classes that can be built without any dependencies
            result = (
                session.exec(
                    (select(AssetClass).join(Recipe).where(Recipe.input_asset_classes == null()))
                )
                .unique()
                .all()
            )
            return result

    def table(self) -> Table:
        """
        Create a table of all asset classes.
        """
        with self._database_session as session:
            asset_classes = session.exec(select(AssetClass)).unique().all()

            # Seek keys are read inside the session context to avoid DetachedInstanceError
            rows = [
                (
                    asset_class.name,
                    asset_class.version,
                    ", ".join(seek_key.name for seek_key in asset_class.seek_keys),
                    asset_class.description,
                )
                for asset_class in asset_classes
            ]

        return build_table(
            "Asset Classes",
            ["Name", "Version", "Seek keys", "Description"],
            rows,
            end_section=True,
        )
