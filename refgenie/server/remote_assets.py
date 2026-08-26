"""
Utility functions for handling remote asset serving.

This module provides common functionality for serving assets from remote locations
when configured, with fallback to local file serving.
"""

import logging
import pathlib

from refgenie.db.tables import Configuration, Remote, RemoteAssetLink
from sqlmodel import Session, select

logger = logging.getLogger(__name__)


def get_remote_url(
    local_path: pathlib.Path,
    session: Session,
    asset_digest: str,
    mode: str,
) -> str | None:
    """
    Find a remote URL for the given local path, if the asset is on a remote.

    Joins through RemoteAssetLink to find remotes that explicitly host this
    (asset_digest, mode) pair. Only returns URLs for assets that have been
    actually pushed (pushed=True). Works for both archive paths (.tgz) and
    individual file paths.

    Args:
        local_path: The local path to the file or archive.
        session: The database session.
        asset_digest: The asset digest to look up.
        mode: The staging mode ("archive" or "file").

    Returns:
        The remote URL as a string if found and constructible, otherwise None.
    """
    query = (
        select(Configuration.genome_stage_folder, Remote.prefix, Remote.type)
        .join(Remote)
        .join(RemoteAssetLink, Remote.id == RemoteAssetLink.remote_id)
        .where(RemoteAssetLink.asset_digest == asset_digest)
        .where(RemoteAssetLink.mode == mode)
        .where(RemoteAssetLink.pushed == True)  # noqa: E712
        .where(Remote.type.in_(["https", "http"]))
        .order_by(Remote.type.desc(), Configuration.id.desc())
    )
    remote_config_data = session.exec(query).first()

    if remote_config_data:
        genome_stage_folder_str, remote_prefix, _ = remote_config_data
        genome_stage_folder = pathlib.Path(genome_stage_folder_str)
        try:
            relative_path = str(local_path.relative_to(genome_stage_folder))
            return f"{remote_prefix.rstrip('/')}/{relative_path.lstrip('/')}"
        except ValueError:
            logger.warning(
                f"Cannot determine remote path for {asset_digest=}. "
                f"{local_path=} is not relative to {genome_stage_folder=}."
            )
    return None
