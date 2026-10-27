"""Data channel handlers for different protocols."""

import ftplib
from pathlib import Path
from typing import Protocol
from urllib.parse import urlparse

from refgenie.logger import logger
from refgenie.utils.http import make_client

__all__ = [
    "ChannelHandler",
    "HTTPChannelHandler",
    "FTPChannelHandler",
    "LocalChannelHandler",
]


class ChannelHandler(Protocol):
    """Protocol for data channel handlers."""

    def fetch_index_content(
        self, index_address: str, credentials: dict | None = None
    ) -> str | None:
        """
        Fetch index.yaml content from the channel.

        Args:
            index_address: Full path/URL to the index.yaml file
            credentials: Optional credentials for authentication

        Returns:
            str | None: Content of index.yaml if successful, None otherwise
        """
        ...

    def test_channel(self, index_address: str, credentials: dict | None = None) -> bool:
        """
        Test channel accessibility.

        Args:
            index_address: Full path/URL to test
            credentials: Optional credentials for authentication

        Returns:
            bool: True if channel is accessible, False otherwise
        """
        ...

    def get_file_url(self, index_address: str, base_dir: str, filename: str) -> str:
        """
        Generate appropriate URL/path for a file.

        Args:
            index_address: Full path/URL to the index.yaml file
            base_dir: Base directory path
            filename: Name of the file

        Returns:
            str: Full URL or path to the file
        """
        ...


class HTTPChannelHandler:
    def fetch_index_content(
        self, index_address: str, credentials: dict | None = None
    ) -> str | None:
        try:
            with make_client(credentials=credentials) as client:
                response = client.get(index_address)
                if response.status_code >= 400:
                    logger.error(f"Failed to download index.yaml from {index_address}")
                    return None
                return response.text

        except Exception as e:
            logger.error(f"Error fetching index.yaml: {str(e)}")
            return None

    def test_channel(self, index_address: str, credentials: dict | None = None) -> bool:
        try:
            with make_client(credentials=credentials) as client:
                response = client.head(index_address)
                return response.status_code < 400
        except Exception as e:
            logger.error(f"HTTP channel test failed: {str(e)}")
            return False

    def get_file_url(self, index_url: str, base_dir: str, filename: str) -> str:
        base_url = index_url.rsplit("/", 1)[0].rstrip("/")
        base_dir = base_dir.strip("/")
        if base_dir:
            base_url = f"{base_url}/{base_dir}"
        return f"{base_url}/{filename}"


class FTPChannelHandler:
    def fetch_index_content(
        self, index_address: str, credentials: dict | None = None
    ) -> str | None:
        try:
            parsed = urlparse(index_address)
            with ftplib.FTP(parsed.netloc) as ftp:
                if credentials:
                    ftp.login(
                        user=credentials.get("username", "anonymous"),
                        passwd=credentials.get("password", "anonymous@"),
                    )
                else:
                    ftp.login()

                content = []
                ftp.retrlines(f"RETR {parsed.path}", content.append)
                return "\n".join(content)

        except Exception as e:
            logger.error(f"Error fetching index.yaml: {str(e)}")
            return None

    def test_channel(self, index_address: str, credentials: dict | None = None) -> bool:
        try:
            parsed = urlparse(index_address)
            with ftplib.FTP(parsed.netloc) as ftp:
                if credentials:
                    ftp.login(
                        user=credentials.get("username", "anonymous"),
                        passwd=credentials.get("password", "anonymous@"),
                    )
                else:
                    ftp.login()
                # Test if index.yaml exists
                ftp.size(parsed.path)  # This will raise error if file doesn't exist
                return True
        except Exception as e:
            logger.error(f"FTP channel test failed: {str(e)}")
            return False

    def get_file_url(self, index_address: str, base_dir: str, filename: str) -> str:
        parsed = urlparse(index_address)
        parent_path = str(Path(parsed.path).parent)
        return f"ftp://{parsed.netloc}{parent_path}/{base_dir}/{filename}"


class LocalChannelHandler:
    def fetch_index_content(
        self, index_address: str, credentials: dict | None = None
    ) -> str | None:
        try:
            path = Path(index_address)
            if not path.exists():
                logger.error(f"index.yaml not found at {path}")
                return None
            return path.read_text()
        except Exception as e:
            logger.error(f"Error fetching index.yaml: {str(e)}")
            return None

    def test_channel(self, index_address: str, credentials: dict | None = None) -> bool:
        try:
            return Path(index_address).exists()
        except Exception as e:
            logger.error(f"Local channel test failed: {str(e)}")
            return False

    def get_file_url(self, index_address: str, base_dir: str, filename: str) -> str:
        parent_path = Path(index_address).parent
        return str(parent_path / base_dir / filename)
