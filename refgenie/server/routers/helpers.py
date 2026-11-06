"""Shared lookups for the server routers: config, staged assets, publication checks."""

import pathlib

from fastapi import HTTPException

from refgenie.core import Refgenie
from refgenie.db.tables import Asset, Configuration, StagedAsset
from refgenie.exceptions import MissingAssetError, MissingConfigDataError
from refgenie.utils.staging import staged_archive_path


def get_config(rgc: Refgenie) -> Configuration | None:
    """Return the server Configuration row in force, or None if absent."""
    try:
        return rgc.configuration.get_latest()
    except MissingConfigDataError:
        return None


def require_config(rgc: Refgenie) -> Configuration:
    """Return the server Configuration, or 500 if it or its genome_stage_folder is missing."""
    config = get_config(rgc)
    if not config or not config.genome_stage_folder:
        raise HTTPException(status_code=500, detail="Server genome_stage_folder not configured")
    return config


def get_asset(rgc: Refgenie, asset_digest: str) -> Asset | None:
    """Return the fully loaded Asset (``serving_modes`` resolvable), or None."""
    try:
        return rgc.asset.get_by_digest(asset_digest)
    except MissingAssetError:
        return None


PUBLISHED_MODES = frozenset({"archive", "file"})


def require_published(
    rgc: Refgenie, asset_digest: str, mode: str, file_path: str | None = None
) -> tuple[Asset, StagedAsset]:
    """Confirm an asset is published in ``mode`` and return ``(asset, staged)``.

    The ``StagedAsset`` row, written only after staging finishes, is the
    publication record; files on disk alone are not trusted. The mode must
    also be one of the asset's current serving modes. A ``file_path`` must be
    listed in the row's ``directory_contents`` (the path-traversal guard).

    Raises:
        HTTPException: 400 for an invalid mode, 404 if any check fails.
    """
    if mode not in PUBLISHED_MODES:
        raise HTTPException(status_code=400, detail=f"Invalid serving mode: {mode}")

    asset = get_asset(rgc, asset_digest)
    if not asset:
        raise HTTPException(status_code=404, detail="Asset not found")

    if mode not in asset.serving_modes:
        raise HTTPException(
            status_code=404, detail=f"Asset {asset_digest} is not served in {mode} mode"
        )

    staged = rgc.stage.find(asset_digest, mode)
    if not staged:
        raise HTTPException(
            status_code=404, detail=f"Asset {asset_digest} is not published in {mode} mode"
        )

    if file_path is not None and (
        mode != "file" or file_path not in (staged.directory_contents or [])
    ):
        raise HTTPException(
            status_code=404, detail=f"File '{file_path}' not found in asset {asset_digest}"
        )

    return asset, staged


def staged_local_path(config: Configuration, asset: Asset, mode: str) -> pathlib.Path:
    """The local path a staged asset is served from.

    Archive mode: the tarball. File mode: the staged directory, a symlink to
    the real files. ``asset.asset_group`` must be loaded.
    """
    if mode == "archive":
        return staged_archive_path(
            config.genome_stage_folder,
            asset.asset_group.genome_digest,
            asset.asset_group.name,
            asset.digest,
        )
    return (
        pathlib.Path(config.genome_stage_folder)
        / asset.asset_group.genome_digest
        / asset.asset_group.name
        / asset.name
    )
