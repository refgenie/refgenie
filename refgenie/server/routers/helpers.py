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

from refgenie.db.tables import Asset, Configuration, StagedAsset
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
        raise HTTPException(
            status_code=500, detail="Server genome_stage_folder not configured"
        )
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
    """Return the Asset with its asset_group eager-loaded, or None."""
    return (
        session.exec(
            select(Asset)
            .options(selectinload(Asset.asset_group))
            .where(Asset.digest == asset_digest)
        )
        .unique()
        .one_or_none()
    )


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
