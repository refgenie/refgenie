from pathlib import Path
from collections.abc import Iterable

from rich.table import Table

from refgenie.db.queries import asset_class as queries
from refgenie.db.queries import recipe as recipe_queries
from refgenie.db.tables import AssetClass, AssetClassSeekKey
from refgenie.exceptions import (
    AssetClassExistsError,
    ConfigError,
    MissingAssetClassError,
)
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.models import Page
from refgenie.utils.io import read_yaml
from refgenie.utils.search import Search
from refgenie.utils.tables import build_table
from refgenie.utils.versioning import validate_semver


class AssetClassManager(ResourceManager):
    """
    A manager for asset classes.
    """

    def add(self, asset_class_source: Path | str, exists_overwrite: bool = False) -> AssetClass:
        """
        Read a asset class from YAML file and register.

        Args:
            asset_class_source: The path to the asset class file or URL.

        Returns:
            AssetClass: The registered asset class.
        """
        logger.debug(f"Loading asset_class from {asset_class_source}")
        return self.add_from_dict(read_yaml(asset_class_source), exists_overwrite=exists_overwrite)

    def add_from_dict(self, asset_class_dict: dict, exists_overwrite: bool = False) -> AssetClass:
        """
        Register an asset class from its definition, in the asset class YAML shape.

        Args:
            asset_class_dict: The definition: ``name``, ``version``,
                ``description``, ``serving_modes`` and ``seek_keys`` (a mapping
                of seek key name to ``value``/``description``/``type``).
            exists_overwrite: Replace an existing class of the same name and
                version instead of raising.

        Returns:
            AssetClass: The registered asset class.
        """
        # Catch a recipe YAML passed by mistake.
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
        with self._database_session as session:
            return queries.get(session, asset_class_name, asset_class_version)

    def get_by_id(self, id: int) -> AssetClass | None:
        """The asset class with this id, with its seek keys, or None if there is none."""
        with self._database_session as session:
            return queries.get_by_id(session, id)

    def page(
        self,
        *,
        name: str | None = None,
        version: str | None = None,
        search: Search | None = None,
        offset: int,
        limit: int,
    ) -> Page[AssetClass]:
        """
        One page of asset classes with their seek keys, by name and version.

        Raises:
            ValueError: If the search names a field that is not searchable.
        """
        with self._database_session as session:
            rows, total = queries.list_page(
                session, name=name, version=version, search=search, offset=offset, limit=limit
            )
        return Page(items=rows, offset=offset, limit=limit, total=total)

    def exists_any_version(self, asset_class_name: str) -> bool:
        """Whether any version of the named asset class is registered."""
        with self._database_session as session:
            return queries.exists_any_version(session, asset_class_name)

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
            return queries.exists(session, asset_class_name, asset_class_version)

    def remove(self, asset_class_name: str, asset_class_version: str | None = None):
        """
        Remove an asset class.

        The class must have no assets and no recipes that use it as input or output.

        If version is omitted and only one version exists, removes it.
        If version is omitted and multiple versions exist, raises an error.

        Args:
            asset_class_name: The name of the asset class.
            asset_class_version: The version of the asset class.
        """
        logger.debug(f"Removing asset class '{asset_class_name}'")
        with self._database_session as session:
            if asset_class_version is not None:
                asset_class = queries.get(
                    session, asset_class_name, asset_class_version, with_groups=True
                )
            else:
                results = queries.list_versions(session, asset_class_name, with_groups=True)
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
            recipe = recipe_queries.find_using_asset_class(session, asset_class.id)
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
            return queries.list_all(session)

    def list_root(self) -> Iterable[AssetClass]:
        """
        List all root asset classes (asset classes that have no dependencies).

        Returns:
            list[AssetClass]: A list of all root asset classes.
        """
        with self._database_session as session:
            return queries.list_root(session)

    def table(self) -> Table:
        """
        Create a table of all asset classes.
        """
        rows = [
            (
                asset_class.name,
                asset_class.version,
                ", ".join(seek_key.name for seek_key in asset_class.seek_keys),
                asset_class.description,
            )
            for asset_class in self.list_all()
        ]

        return build_table(
            "Asset Classes",
            ["Name", "Version", "Seek keys", "Description"],
            rows,
            end_section=True,
        )
