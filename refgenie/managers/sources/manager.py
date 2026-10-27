"""
SourceManager - data channels, where recipes and asset classes come from.

Refgenieservers that assets are pulled from are ``ServerManager``'s job
(``servers.py``), not this one's.
"""

from pathlib import Path
from collections.abc import Iterator
from typing import TYPE_CHECKING, Literal

import yaml
from pydantic import BaseModel, HttpUrl
from rich.table import Table
from sqlalchemy.engine import Engine
from sqlmodel import select

from refgenie.db.tables import Configuration, DataChannel, DataChannelType
from refgenie.exceptions import AssetClassExistsError, RecipeExistsError, RefgenieError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.sources.handlers import (
    ChannelHandler,
    FTPChannelHandler,
    HTTPChannelHandler,
    LocalChannelHandler,
)
from refgenie.utils.encryption import decrypt_credentials, encrypt_dict
from refgenie.models import DataChannelSyncReport
from refgenie.utils.tables import build_table

if TYPE_CHECKING:
    from refgenie.managers.asset_class import AssetClassManager
    from refgenie.managers.recipe import RecipeManager


class IndexFileSection(BaseModel):
    dir: str = ""
    files: list[HttpUrl | Path]


class IndexFile(BaseModel):
    """
    Pydantic schema for validating index files fetched from data channels.
    """

    asset_class: IndexFileSection
    recipe: IndexFileSection


class SourceManager(ResourceManager):
    """Data channels: registering them, reading the recipes and asset classes they
    list, and syncing those into the local asset-class and recipe managers."""

    def __init__(
        self,
        database_engine: Engine,
        asset_class_manager: "AssetClassManager",
        recipe_manager: "RecipeManager",
    ):
        """
        Initialize the SourceManager.

        Args:
            database_engine: The database engine.
            asset_class_manager: Where ``sync_channel`` registers asset classes.
            recipe_manager: Where ``sync_channel`` registers recipes.
        """
        super().__init__(database_engine)
        self._asset_class = asset_class_manager
        self._recipe = recipe_manager
        self._handlers: dict[DataChannelType, ChannelHandler] = {
            DataChannelType.http: HTTPChannelHandler(),
            DataChannelType.https: HTTPChannelHandler(),
            DataChannelType.ftp: FTPChannelHandler(),
            DataChannelType.local: LocalChannelHandler(),
        }

    # =========================================================================
    # Data Channel CRUD Operations
    # =========================================================================

    def _get_handler(self, channel_type: DataChannelType) -> ChannelHandler:
        """
        Get appropriate handler for the given channel type.

        Args:
            channel_type: Type of the data channel

        Returns:
            ChannelHandler: Handler for the given channel type
        """
        return self._handlers[channel_type]

    def add_channel(
        self,
        name: str,
        type: DataChannelType,
        index_address: str,
        description: str | None = None,
        credentials: dict | None = None,
    ) -> DataChannel:
        """
        Add a new data channel to the database.

        Args:
            name: Unique name for the channel
            type: Type of the channel (FTP, HTTP, etc.)
            index_address: Full path/URL to the index.yaml file
            description: Optional description of the channel
            credentials: Optional credentials for authenticated access

        Returns:
            DataChannel: The created data channel object

        Raises:
            ValueError: If a data channel with the given name already exists
        """
        # double underscores are not allowed in the name as they may be used
        # to separate the channel from the file name
        if "__" in name:
            raise ValueError("Channel name cannot contain '__'")

        with self._database_session as session:
            # Queried through this session: get_channel opens its own, and
            # _database_session must not be re-entered from inside an open block.
            if session.exec(select(DataChannel).where(DataChannel.name == name)).first():
                raise ValueError(f"Data channel '{name}' already exists")

            configuration = session.exec(select(Configuration)).unique().one()
            channel = DataChannel(
                name=name,
                type=type,
                index_address=index_address,
                description=description,
                encrypted_credentials=encrypt_dict(credentials),
                configuration_id=configuration.id,
            )
            session.add(channel)
            session.commit()
            session.refresh(channel)
            return channel

    def remove_channel(self, name: str) -> bool:
        """
        Remove a data channel from the database.

        Args:
            name: Name of the channel to remove

        Returns:
            bool: True if channel was removed, False if not found
        """
        with self._database_session as session:
            channel = session.exec(select(DataChannel).where(DataChannel.name == name)).first()
            if channel:
                session.delete(channel)
                session.commit()
                return True
            return False

    def get_channel(self, name: str) -> DataChannel | None:
        """
        Get a data channel by name.

        Args:
            name: Name of the channel to retrieve

        Returns:
            DataChannel | None: The channel if found, None otherwise
        """
        with self._database_session as session:
            return session.exec(select(DataChannel).where(DataChannel.name == name)).first()

    def list_channels(self) -> list[DataChannel]:
        """
        List all configured data channels.

        Returns:
            list[DataChannel]: List of all data channels
        """
        with self._database_session as session:
            return session.exec(select(DataChannel)).unique().all()

    def test_channel(self, name: str) -> bool:
        """
        Test if a data channel is accessible.

        Args:
            name: Name of the channel to test

        Returns:
            bool: True if channel is accessible, False otherwise
        """
        if not (channel := self.get_channel(name)):
            return False

        handler = self._get_handler(channel.type)

        if not handler.test_channel(
            channel.index_address, decrypt_credentials(channel.encrypted_credentials)
        ):
            logger.error(f"Failed to access channel '{name}'")
            return False

        try:
            assert self.get_index_file(name) is not None, "Index file is empty"
        except Exception as e:
            logger.error(f"Failed to fetch index file for channel '{name}': {e}")
            return False
        return True

    def get_index_file(self, name: str) -> IndexFile | None:
        """
        Get the index file content for a data channel.

        Args:
            name: Name of the channel to retrieve index file from

        Returns:
            IndexFile | None: Index file content if found, None otherwise
        """
        if not (channel := self.get_channel(name)):
            logger.error(f"Channel '{name}' not found")
            return None

        handler = self._get_handler(channel.type)

        if not (
            content := handler.fetch_index_content(
                channel.index_address,
                decrypt_credentials(channel.encrypted_credentials),
            )
        ):
            return None

        parsed = yaml.safe_load(content)
        if not isinstance(parsed, dict):
            logger.error(
                f"Index file for channel '{name}' is not valid YAML mapping "
                f"(got {type(parsed).__name__}). "
                f"Check that the channel's index_address points to an index.yaml file."
            )
            return None

        return IndexFile(**parsed)

    def _iter_files(
        self,
        channel_name: str,
        section: Literal["asset_class", "recipe"],
    ) -> Iterator[str]:
        """
        Generic iterator for files listed in index file of the channel.

        Args:
            channel_name: Name of the channel containing the index file
            section: Section name in index ('asset_class' or 'recipe')

        Yields:
            str: URLs or paths to the files listed in the index
        """
        if not (channel := self.get_channel(channel_name)):
            logger.error(f"Channel '{channel_name}' not found")
            return

        if (index_file := self.get_index_file(channel_name)) is None:
            logger.error(f"Failed to fetch index file for channel '{channel_name}'")
            return

        files = (
            index_file.asset_class.files if section == "asset_class" else index_file.recipe.files
        )
        files_dir = (
            index_file.asset_class.dir if section == "asset_class" else index_file.recipe.dir
        )

        handler = self._get_handler(channel.type)

        for file in files:
            if isinstance(file, HttpUrl):
                yield str(file)
            else:
                yield handler.get_file_url(channel.index_address, files_dir, file)

    def iter_asset_classes(self, channel_name: str) -> Iterator[str]:
        """
        Iterator for asset class files listed in index file

        Args:
            channel_name: Name of the channel containing the indexfile

        Yields:
            str: URLs or paths to the asset class files
        """
        yield from self._iter_files(channel_name, "asset_class")

    def iter_recipes(self, channel_name: str) -> Iterator[str]:
        """
        Iterator for recipe files listed in index file

        Args:
            channel_name: Name of the channel containing the index file

        Yields:
            str: URLs or paths to the recipe files
        """
        yield from self._iter_files(channel_name, "recipe")

    def channels_table(self, channel_names: list[str] | None = None) -> Table:
        """
        Create a formatted table of all data channels

        Args:
            channel_names: Optional list of channel names to include in the table

        Returns:
            Table: Rich table object with data channel information
        """
        return build_table(
            "Data Channels",
            [("Name", "cyan"), "Type", "Index Address", "Description", "Credentials set"],
            [
                (
                    channel.name,
                    channel.type.value,
                    channel.index_address,
                    channel.description or "",
                    str(bool(channel.encrypted_credentials)),
                )
                for channel in self.list_channels()
                if not channel_names or channel.name in channel_names
            ],
        )

    def sync_channel(
        self,
        channel_name: str,
        *,
        exists_ok: bool = False,
        exists_overwrite: bool = False,
    ) -> DataChannelSyncReport:
        """
        Register every asset class and recipe a data channel publishes.

        The channel is read here; the items land in the injected asset-class
        and recipe managers. Asset classes go first, and a failure there stops
        the recipes: a recipe whose asset class did not register would only
        fail again, louder.

        Args:
            channel_name: The data channel to sync from.
            exists_ok: Skip items already registered instead of failing them.
            exists_overwrite: Replace items already registered.

        Returns:
            DataChannelSyncReport: Counts and per-item errors. ``report.ok`` is
            False when anything failed; nothing is raised for item failures.

        Raises:
            RefgenieError: If the channel does not exist or cannot be reached.
        """
        if not self.test_channel(channel_name):
            raise RefgenieError(f"Data channel '{channel_name}' is not accessible")

        report = DataChannelSyncReport(channel=channel_name)
        for asset_class_url in self.iter_asset_classes(channel_name):
            try:
                self._asset_class.add(asset_class_url, exists_overwrite=exists_overwrite)
                report.asset_classes_added += 1
            except AssetClassExistsError as e:
                if exists_ok:
                    report.asset_classes_skipped += 1
                    continue
                report.asset_classes_failed += 1
                report.errors.append(f"asset class {asset_class_url}: {e}")
            except (RefgenieError, OSError, ValueError) as e:
                report.asset_classes_failed += 1
                report.errors.append(f"asset class {asset_class_url}: {e}")

        if report.asset_classes_failed:
            return report

        for recipe_url in self.iter_recipes(channel_name):
            try:
                self._recipe.add(recipe_url, exists_overwrite=exists_overwrite)
                report.recipes_added += 1
            except RecipeExistsError as e:
                if exists_ok:
                    report.recipes_skipped += 1
                    continue
                report.recipes_failed += 1
                report.errors.append(f"recipe {recipe_url}: {e}")
            except (RefgenieError, OSError, ValueError) as e:
                report.recipes_failed += 1
                report.errors.append(f"recipe {recipe_url}: {e}")
        return report
