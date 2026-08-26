"""
SourceManager - manages all external data sources.

This combines:
- Refgenieserver connections (for pulling assets, resolving aliases remotely)
- Data channels (for fetching recipes and asset classes from external sources)
"""

from pathlib import Path
from collections.abc import Iterator
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from refgenie.managers.configuration import ConfigurationManager
import yaml
from pydantic import BaseModel, HttpUrl
from rich.table import Table
from sqlalchemy.engine import Engine
from sqlmodel import select

from refgenie.db.tables import Configuration, DataChannel, DataChannelType
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.sources.handlers import (
    ChannelHandler,
    FTPChannelHandler,
    HTTPChannelHandler,
    LocalChannelHandler,
)
from refgenie.managers.sources.client import RefgenieserverClient, ServerClient
from refgenie.utils.encryption import decrypt_credentials, encrypt_dict
from refgenie.utils.tables import build_table


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
    """
    Manages all external data sources:
    - Refgenieserver connections (for assets, genomes, aliases)
    - Data channels (for recipes, asset classes)
    """

    def __init__(
        self,
        database_engine: Engine,
        configuration_manager: "ConfigurationManager",  # type: ignore  # Forward reference
        server_clients: dict[str, ServerClient | None] = None,
    ):
        """
        Initialize the SourceManager.

        Args:
            database_engine: The database engine.
            configuration_manager: The ConfigurationManager for server subscriptions.
            server_clients: Optional pre-configured server clients mapping.
        """
        super().__init__(database_engine)
        self._config = configuration_manager
        self._server_clients: dict[str, ServerClient] = server_clients or {}
        self._handlers: dict[DataChannelType, ChannelHandler] = {
            DataChannelType.http: HTTPChannelHandler(),
            DataChannelType.https: HTTPChannelHandler(),
            DataChannelType.ftp: FTPChannelHandler(),
            DataChannelType.local: LocalChannelHandler(),
        }

    # =========================================================================
    # Server Client Operations
    # =========================================================================

    @property
    def server_clients(self) -> dict[str, ServerClient]:
        """Get the server clients mapping."""
        return self._server_clients

    def get_server_client(self, url: str) -> ServerClient:
        """
        Get or create a client for a refgenieserver URL.

        Args:
            url: The server URL.

        Returns:
            ServerClient: The server client.
        """
        if url not in self._server_clients:
            self._server_clients[url] = RefgenieserverClient(url)
        return self._server_clients[url]

    def get_subscriptions(self) -> list[str]:
        """Get the list of subscribed server URLs."""
        return list(self._config.get_server_subscriptions())

    def sync_server_clients(
        self, server_clients_override: dict[str, ServerClient | None] = None
    ) -> None:
        """
        Sync server clients with current subscriptions, or override with provided mapping.

        Args:
            server_clients_override: The server clients mapping to override.
        """
        if server_clients_override is not None:
            self._server_clients = server_clients_override
        else:
            for url in self._config.get_server_subscriptions():
                if url not in self._server_clients:
                    self._server_clients[url] = RefgenieserverClient(url)

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
