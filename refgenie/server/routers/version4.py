"""Version 4 API router: archive and file downloads plus summary statistics.

Server mode only; listings shared with local mode are in ``catalog``.
"""

import logging
import pathlib

from fastapi import APIRouter, Depends, HTTPException, Query

from refgenie.core import Refgenie
from starlette.responses import FileResponse, RedirectResponse

from refgenie.server.dependencies import get_refgenie
from refgenie.server.routers.helpers import (
    require_config,
    require_published,
    staged_local_path,
)
from refgenie.server.schemas import (
    ArchiveRecord,
    DatabaseSummaryResponse,
    SpeciesStatistics,
)
from refgenie.server.pagination import PageParams, PageResponse, get_page_params, respond

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/archives", response_model=PageResponse[ArchiveRecord])
def list_archives(
    digest: str | None = Query(None, description="The asset digest"),
    asset_digest: str | None = Query(None, description="The asset digest"),
    genome_digest: str | None = Query(
        None, description="Return only archives belonging to this genome"
    ),
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """List downloadable archives (staged assets in archive mode).

    ``digest`` and ``asset_digest`` both filter by asset digest. Filters
    combine with AND. Only archives the download route will serve are listed.
    """
    page = rgc.stage.archives(
        asset_digest=asset_digest or digest,
        genome_digest=genome_digest,
        offset=paging.offset,
        limit=paging.page_size,
    )
    return respond(
        page,
        paging,
        lambda sa: ArchiveRecord(
            digest=sa.asset_digest,
            asset_digest=sa.asset_digest,
            tarball_digest=sa.tarball_digest,
            size=sa.tarball_size,
            directory_contents=sa.directory_contents,
            build_commands=sa.build_commands,
            download_count=sa.download_count,
        ),
    )


@router.get("/archives/{asset_digest}/download", response_class=FileResponse)
def download_archive(
    asset_digest: str,
    rgc: Refgenie = Depends(get_refgenie),
):
    """Download an archive by its asset digest.

    Redirects to a remote the archive was pushed to, if any; otherwise serves
    the local file.
    """
    asset, _ = require_published(rgc, asset_digest, "archive")

    config = require_config(rgc)
    tarball_path = staged_local_path(config, asset, "archive")

    if urls := rgc.remote.pushed_urls(asset_digest, {"archive": tarball_path}):
        remote_url = urls[0][2]
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
    rgc: Refgenie = Depends(get_refgenie),
):
    """Download one file of an asset.

    The file must be listed in the staged ``directory_contents``. Redirects to
    a remote the file was pushed to, if any; otherwise serves it from
    genome_stage_folder.
    """
    asset, _ = require_published(rgc, asset_digest, "file", file_path)

    config = require_config(rgc)

    local_file_path = staged_local_path(config, asset, "file") / file_path

    if urls := rgc.remote.pushed_urls(asset_digest, {"file": local_file_path}):
        remote_url = urls[0][2]
        logger.debug(f"Redirecting file request to remote URL: {remote_url}")
        return RedirectResponse(url=remote_url)

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
    """Healthcheck for monitoring."""
    return {"status": "ok"}


@router.get("/species/summary", response_model=dict[str, SpeciesStatistics])
def species_summary(rgc: Refgenie = Depends(get_refgenie)) -> dict[str, SpeciesStatistics]:
    """Counts of genomes, asset classes and assets per species.

    Only genomes with assets are counted. Genomes with no species name are
    grouped under "Unknown".
    """
    return {
        species: SpeciesStatistics(**counts)
        for species, counts in rgc.database.species_summary().items()
    }


@router.get("/summary", response_model=DatabaseSummaryResponse)
def summary(rgc: Refgenie = Depends(get_refgenie)) -> DatabaseSummaryResponse:
    """Total counts of genomes, asset groups and assets."""
    return DatabaseSummaryResponse(**rgc.database.summary())
