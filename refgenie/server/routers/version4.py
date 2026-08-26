"""
Version 4 API router for refgenieserver.

Handlers take their dependencies (database session) through FastAPI's
``Depends`` from ``refgenie.server.dependencies``, so ``app.dependency_overrides``
works and there is no module-level state. This router is included in server mode
only; entity listings shared with local mode live in
``refgenie.server.routers.shared``. This router owns the archive and
file-serving endpoints plus the summary statistics.
"""

import logging
import pathlib

from fastapi import APIRouter, Depends, HTTPException, Query
from refgenie.db.tables import (
    StagedAsset,
    Asset,
    AssetClass,
    AssetGroup,
    Genome,
)
from sqlmodel import Session, func, select
from starlette.responses import FileResponse, RedirectResponse

from refgenie.server.dependencies import get_db_session
from refgenie.server.routers.helpers import (
    get_asset_with_group,
    get_staged_asset,
    require_config,
    staged_local_path,
)
from refgenie.server.schemas import (
    ArchiveRecord,
    DatabaseSummaryResponse,
    SpeciesStatistics,
)
from refgenie.server.remote_assets import get_remote_url
from refgenie.utils.pagination import (
    PaginatedResponse,
    PaginationParams,
    get_pagination_params,
    paginate_list,
)

# Configure logger
logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/archives", response_model=PaginatedResponse[ArchiveRecord])
def list_archives(
    digest: str | None = Query(None, description="The asset digest"),
    asset_digest: str | None = Query(None, description="The asset digest"),
    genome_digest: str | None = Query(
        None, description="Return only archives belonging to this genome"
    ),
    pagination: PaginationParams = Depends(get_pagination_params),
    session: Session = Depends(get_db_session),
):
    """
    List downloadable archives (staged assets in archive mode).

    ``digest`` and ``asset_digest`` are interchangeable filters -- both match
    the asset digest, which is the identifier the download route uses.
    ``genome_digest`` restricts results to archives of a single genome (joins
    StagedAsset -> Asset -> AssetGroup). All filters combine with AND.
    """
    query = select(StagedAsset).where(StagedAsset.mode == "archive")
    wanted = asset_digest or digest
    if wanted:
        query = query.where(StagedAsset.asset_digest == wanted)
    if genome_digest:
        query = (
            query.join(Asset, Asset.digest == StagedAsset.asset_digest)
            .join(AssetGroup, AssetGroup.id == Asset.asset_group_id)
            .where(AssetGroup.genome_digest == genome_digest)
        )
    staged = session.exec(query).all()
    records = [
        ArchiveRecord(
            digest=sa.asset_digest,
            asset_digest=sa.asset_digest,
            tarball_digest=sa.tarball_digest,
            size=sa.tarball_size,
            directory_contents=sa.directory_contents,
            build_commands=sa.build_commands,
            download_count=sa.download_count,
        )
        for sa in staged
    ]
    return paginate_list(records, pagination)


@router.get("/archives/{asset_digest}/download", response_class=FileResponse)
def download_archive(asset_digest: str, session: Session = Depends(get_db_session)):
    """
    Download an archive by its asset digest.

    The archive can be served in two ways:
    1.  **Remote Redirect**: If a suitable remote configuration (type "http" or "https")
        is found and the archive path is relative to the server's genome archive folder,
        the endpoint redirects to the archive's URL.
    2.  **Local File Serve**: If no suitable remote configuration is found, the archive
        is served directly from the local filesystem.
    """
    sa = get_staged_asset(session, asset_digest, "archive")

    if not sa:
        raise HTTPException(
            status_code=404,
            detail=f"No archive available for asset {asset_digest}. "
            f"Use /assets/{asset_digest}/files/ for file-level access.",
        )

    # Derive tarball path from convention
    asset = get_asset_with_group(session, asset_digest)

    if not asset:
        raise HTTPException(status_code=404, detail=f"Asset {asset_digest} not found")

    config = require_config(session)
    tarball_path = staged_local_path(config, asset, "archive")

    if remote_url := get_remote_url(tarball_path, session, asset_digest, "archive"):
        logger.debug(f"Redirecting to remote URL: {remote_url}")
        return RedirectResponse(url=remote_url)

    logger.debug(f"Serving archive {asset_digest} from local path: {tarball_path}")
    if not tarball_path.is_file():
        logger.exception(f"Archive file not found at path: {tarball_path}")
        raise HTTPException(status_code=404, detail=f"Archive file not found for {asset_digest}")

    return FileResponse(
        tarball_path,
        filename=tarball_path.name,
        media_type="application/octet-stream",
    )


@router.get("/assets/{asset_digest}/files/{file_path:path}", response_class=FileResponse)
def download_asset_file(
    asset_digest: str,
    file_path: str,
    session: Session = Depends(get_db_session),
):
    """
    Download an individual file from an asset.

    The file must be listed in the StagedAsset's directory_contents (prevents path traversal).
    Files are served from genome_stage_folder -- for file-mode assets, this is a symlink
    that the OS follows transparently to the real file in genome_folder.
    Supports remote redirect via get_remote_url() if configured.
    """
    # Verify the asset is staged for file serving
    sa = get_staged_asset(session, asset_digest, "file")

    if not sa:
        raise HTTPException(
            status_code=404, detail=f"Asset {asset_digest} is not staged for file-level serving"
        )

    # Path traversal protection: file_path must be in directory_contents
    if file_path not in sa.directory_contents:
        raise HTTPException(
            status_code=404, detail=f"File '{file_path}' not found in asset {asset_digest}"
        )

    # Derive local path from convention
    asset = get_asset_with_group(session, asset_digest)

    if not asset:
        raise HTTPException(status_code=404, detail=f"Asset {asset_digest} not found")

    config = require_config(session)

    # Resolve from genome_stage_folder -- for file-mode, this path is a symlink
    # that the OS follows to the real file in genome_folder
    local_file_path = staged_local_path(config, asset, "file") / file_path

    # Check for remote URL redirect
    if remote_url := get_remote_url(local_file_path, session, asset_digest, "file"):
        logger.debug(f"Redirecting file request to remote URL: {remote_url}")
        return RedirectResponse(url=remote_url)

    # Serve locally
    if not local_file_path.is_file():
        logger.error(f"File not found at path: {local_file_path}")
        raise HTTPException(status_code=404, detail=f"File not found on disk: {file_path}")

    return FileResponse(
        local_file_path,
        filename=pathlib.Path(file_path).name,
        media_type="application/octet-stream",
    )


@router.get("/healthcheck", include_in_schema=False)
def healthcheck():
    """
    Simple healthcheck endpoint for monitoring.
    Returns 200 OK and a JSON status message.
    """
    return {"status": "ok"}


@router.get("/species/summary", response_model=dict[str, SpeciesStatistics])
def species_summary(
    session: Session = Depends(get_db_session),
) -> dict[str, SpeciesStatistics]:
    """
    Returns a summary of genomes, asset classes, and assets grouped by species.

    This endpoint aggregates data for each species showing:
    - Number of genomes
    - Number of unique asset classes
    - Total number of assets

    Returns:
        dict[str, SpeciesStatistics]: Dictionary with species names as keys and
        SpeciesStatistics objects as values. Species with null names are grouped
        under "Unknown".
    """
    # Query to get counts grouped by species
    result = session.exec(
        select(
            Genome.species_name,
            func.count(func.distinct(Genome.digest)).label("genome_count"),
            func.count(func.distinct(AssetClass.id)).label("asset_class_count"),
            func.count(Asset.digest).label("asset_count"),
        )
        .select_from(Genome)
        .join(AssetGroup, Genome.digest == AssetGroup.genome_digest)
        .join(AssetClass, AssetGroup.asset_class_id == AssetClass.id)
        .join(Asset, AssetGroup.id == Asset.asset_group_id)
        .group_by(Genome.species_name)
    ).all()

    # Format the response
    summary = {}
    for row in result:
        species = row.species_name if row.species_name else "Unknown"
        summary[species] = SpeciesStatistics(
            genomes=row.genome_count,
            asset_classes=row.asset_class_count,
            assets=row.asset_count,
        )

    return summary


@router.get("/summary", response_model=DatabaseSummaryResponse)
def summary(
    session: Session = Depends(get_db_session),
) -> DatabaseSummaryResponse:
    """
    Returns a summary of the Refgenieserver database.

    This endpoint provides counts of:
    - Total genomes
    - Total asset groups
    - Total assets

    Returns:
        DatabaseSummaryResponse: Object with summary counts.
    """
    total_genomes = session.exec(select(func.count(func.distinct(Genome.digest)))).one()
    total_asset_groups = session.exec(select(func.count(func.distinct(AssetGroup.id)))).one()
    total_assets = session.exec(select(func.count(func.distinct(Asset.digest)))).one()

    return DatabaseSummaryResponse(
        genomes=total_genomes,
        asset_groups=total_asset_groups,
        assets=total_assets,
    )
