"""
Shared data-access helpers for the server routers.

These consolidate the lookup patterns repeated across ``version4``,
``ga4gh_drs`` and ``shared``: the server Configuration row (with its
``genome_stage_folder`` guard), StagedAsset records keyed by asset digest and
serving mode, the Asset-with-asset_group eager load, and the convention that
maps a staged asset to its on-disk location (tarball for archive mode,
symlinked directory for file mode).
"""

import pathlib

from fastapi import HTTPException
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.tables import Asset, AssetGroup, Configuration, StagedAsset
from refgenie.utils.staging import staged_archive_path


def get_config(session: Session) -> Configuration | None:
    """Return the server Configuration row, or None if absent."""
    return session.exec(select(Configuration)).first()


def require_config(session: Session) -> Configuration:
    """
    Return the server Configuration, raising 500 if it (or its
    genome_stage_folder) is missing -- serving endpoints cannot derive local
    paths without a stage folder.
    """
    config = get_config(session)
    if not config or not config.genome_stage_folder:
        raise HTTPException(status_code=500, detail="Server genome_stage_folder not configured")
    return config


def get_staged_asset(session: Session, asset_digest: str, mode: str) -> StagedAsset | None:
    """Return the StagedAsset for (asset_digest, mode), or None if not staged."""
    return session.exec(
        select(StagedAsset).where(
            StagedAsset.asset_digest == asset_digest,
            StagedAsset.mode == mode,
        )
    ).one_or_none()


def get_asset_with_group(session: Session, asset_digest: str) -> Asset | None:
    """
    Return the Asset with its asset_group -> asset_class eager-loaded, or None.

    The asset_class load is what lets ``asset.serving_modes`` resolve without
    a lazy load (it falls back to the asset class's serving_modes when the
    asset has no override).
    """
    return (
        session.exec(
            select(Asset)
            .options(selectinload(Asset.asset_group).selectinload(AssetGroup.asset_class))
            .where(Asset.digest == asset_digest)
        )
        .unique()
        .one_or_none()
    )


PUBLISHED_MODES = frozenset({"archive", "file"})


def require_published(
    session: Session, asset_digest: str, mode: str, file_path: str | None = None
) -> tuple[Asset, StagedAsset]:
    """
    Confirm an asset is published in ``mode`` and return ``(asset, staged)``.

    A ``StagedAsset`` row is the publication record: it is written only after
    staging finishes writing (and, for archive mode, checksumming) the bytes
    on disk. A file sitting at the predictable stage path with no matching row
    is never evidence of publication -- it may be a leftover from an
    interrupted or superseded stage -- so every byte-serving and access-URL
    route must go through this check rather than trusting the filesystem.

    Checks, in order:
    - ``mode`` must be one of :data:`PUBLISHED_MODES` -- 400 otherwise.
    - The asset must exist -- 404 otherwise.
    - ``mode`` must be one of the asset's effective serving modes -- 404
      otherwise (an asset class or override change can un-publish a mode even
      though an old StagedAsset row and tarball still exist).
    - A ``StagedAsset(asset_digest, mode)`` row must exist -- 404 otherwise.
    - If ``file_path`` is given, ``mode`` must be "file" and ``file_path``
      must be in the staged row's ``directory_contents`` -- 404 otherwise
      (this is the path-traversal guard for file-level access).

    Raises:
        HTTPException: 400 for an invalid mode, 404 for anything else above.
    """
    if mode not in PUBLISHED_MODES:
        raise HTTPException(status_code=400, detail=f"Invalid serving mode: {mode}")

    asset = get_asset_with_group(session, asset_digest)
    if not asset:
        raise HTTPException(status_code=404, detail="Asset not found")

    if mode not in asset.serving_modes:
        raise HTTPException(
            status_code=404, detail=f"Asset {asset_digest} is not served in {mode} mode"
        )

    staged = get_staged_asset(session, asset_digest, mode)
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
    """
    Derive the local path a staged asset is served from.

    Archive mode: the tarball path (content-addressed, via staged_archive_path).
    File mode: the staged asset directory -- a symlink the OS follows to the
    real files in genome_folder; callers append individual file names.

    ``asset.asset_group`` must be loaded (use :func:`get_asset_with_group`).
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
