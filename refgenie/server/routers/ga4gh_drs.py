import datetime
import logging

from fastapi import APIRouter, Depends, HTTPException, Request
from fastapi.responses import FileResponse
from starlette.responses import RedirectResponse

from refgenie.core import Refgenie
from refgenie.db.tables import Asset
from refgenie.server.const import (
    SERVICE_CONTACT_URL,
    SERVICE_DOCUMENTATION_URL,
    SERVICE_ORGANIZATION,
)
from refgenie.server.dependencies import get_refgenie
from refgenie.server.routers.helpers import (
    PUBLISHED_MODES,
    get_asset,
    get_config,
    require_config,
    require_published,
    staged_local_path,
)
from refgenie.server.schemas import (
    AccessMethod,
    AccessURL,
    Checksum,
    ContentsObject,
    DRSModel,
    Organization,
    ServiceInfo,
    ServiceType,
)

router = APIRouter()

DRS_SPEC_VERSION = "1.2.0"
logger = logging.getLogger(__name__)


def parse_refgenie_drs_object_id(object_id: str) -> tuple[str, str | None]:
    """Parse a DRS object ID, ``{digest}`` or ``{digest}:{filename}``, into (digest, file_path)."""
    if not object_id or len(object_id) < 8:
        raise HTTPException(status_code=400, detail=f"Object ID {object_id} is malformed")

    if ":" in object_id:
        digest, file_path = object_id.split(":", 1)
        return digest, file_path
    return object_id, None


def _archive_checksums(rgc: Refgenie, asset: Asset) -> list[Checksum]:
    """DRS ``checksums`` for the downloaded bytes: the archive tarball's sha-256.

    Not ``asset.digest``, which identifies content but is not a byte checksum.
    Empty when there is no staged archive.
    """
    staged_archive = rgc.stage.find(asset.digest, "archive")
    if staged_archive and staged_archive.tarball_digest:
        return [Checksum(checksum=staged_archive.tarball_digest, type="sha-256")]
    return []


def _build_access_methods(
    rgc: Refgenie,
    asset: Asset,
    base_uri: str,
    object_id: str,
    modes: set[str],
    file_path: str | None = None,
) -> list[AccessMethod]:
    """DRS access methods for every published, served mode, local and remote.

    access_id scheme:
    - "archive"                              -- local archive (this server)
    - "file"                                 -- local files (this server)
    - "{remote_type}:{remote_id}:archive"    -- remote archive (e.g., "s3:1:archive")
    - "{remote_type}:{remote_id}:file"       -- remote files (e.g., "https:2:file")
    """
    methods = []

    # Only modes that are both served and staged; a serving-mode change can
    # leave stale staged copies behind.
    published = {m for m in modes & PUBLISHED_MODES if rgc.stage.find(asset.digest, m)}

    config = get_config(rgc)
    has_local = config and config.genome_stage_folder

    if has_local:
        if "archive" in published:
            methods.append(
                AccessMethod(
                    type="https",
                    access_url=AccessURL(
                        url=f"{base_uri}/ga4gh/drs/v1/objects/{object_id}/access/archive",
                        headers=None,
                    ),
                    access_id="archive",
                )
            )

        if "file" in published:
            methods.append(
                AccessMethod(
                    type="https",
                    access_url=AccessURL(
                        url=f"{base_uri}/ga4gh/drs/v1/objects/{object_id}/access/file",
                        headers=None,
                    ),
                    access_id="file",
                )
            )

    # Remote URLs are the staged path re-rooted at the remote, so they need the stage folder too.
    if has_local and published:
        asset_with_group = get_asset(rgc, asset.digest)
        if asset_with_group:
            local_paths = {
                mode: staged_local_path(config, asset_with_group, mode) for mode in published
            }
            if file_path is not None and "file" in local_paths:
                local_paths["file"] = local_paths["file"] / file_path
            for remote, mode, remote_url in rgc.remote.pushed_urls(asset.digest, local_paths):
                methods.append(
                    AccessMethod(
                        type=remote.type.value,
                        access_url=AccessURL(url=remote_url, headers=None),
                        access_id=f"{remote.type.value}:{remote.id}:{mode}",
                    )
                )

    return methods


@router.get("/service-info", response_model=ServiceInfo)
async def get_service_info(request: Request):
    """Return information about the DRS service."""
    return ServiceInfo(
        id="org.refgenie.server",
        name="Refgenie Server DRS",
        type=ServiceType(group="org.ga4gh", artifact="drs", version=DRS_SPEC_VERSION),
        description=(
            "Refgenie Server implementing GA4GH DRS specification for genomic reference assets"
        ),
        organization=Organization(**SERVICE_ORGANIZATION),
        contactUrl=SERVICE_CONTACT_URL,
        documentationUrl=SERVICE_DOCUMENTATION_URL,
        createdAt=datetime.datetime.now(),
        updatedAt=datetime.datetime.now(),
        environment="production",
        version=DRS_SPEC_VERSION,
    )


@router.get("/objects/{object_id}", response_model=DRSModel)
async def get_drs_object_metadata(
    object_id: str,
    request: Request,
    rgc: Refgenie = Depends(get_refgenie),
):
    """Return metadata about a DrsObject.

    File mode gives a bundle, archive-only a simple object, and no mode
    metadata only. A ``{digest}:{filename}`` id gives a single-file object.
    """
    asset_digest, file_path = parse_refgenie_drs_object_id(object_id)
    base_uri = f"{request.url.scheme}://{request.url.netloc}"

    asset = get_asset(rgc, asset_digest)

    if not asset:
        raise HTTPException(status_code=404, detail=f"Asset with digest {asset_digest} not found")

    modes = set(asset.serving_modes)

    # --- File-level request ({digest}:{filename}) ---
    if file_path is not None:
        require_published(rgc, asset_digest, "file", file_path)

        access_methods = _build_access_methods(
            rgc, asset, base_uri, object_id, modes={"file"}, file_path=file_path
        )

        return DRSModel(
            id=object_id,
            name=file_path,
            self_uri=f"drs://{request.url.netloc}/{object_id}",
            size=0,  # per-file sizes are not tracked
            created_time=asset.created_at or datetime.datetime.now(),
            updated_time=asset.created_at or datetime.datetime.now(),
            checksums=[],  # no per-file byte digest is tracked
            access_methods=access_methods,
            contents=None,
            description=f"File: {file_path}",
            aliases=[],
        )

    # File mode (with or without archive): a DRS bundle
    if "file" in modes:
        staged_file = rgc.stage.find(asset.digest, "file")

        contents = [
            ContentsObject(
                name=filename,
                id=f"{asset.digest}:{filename}",
                drs_uri=[f"drs://{request.url.netloc}/{asset.digest}:{filename}"],
                contents=None,
            )
            for filename in (staged_file.directory_contents if staged_file else [])
        ]

        access_methods = _build_access_methods(rgc, asset, base_uri, object_id, modes)

        return DRSModel(
            id=object_id,
            name=asset.registry_path,
            self_uri=f"drs://{request.url.netloc}/{object_id}",
            size=asset.size or 0,
            created_time=asset.created_at or datetime.datetime.now(),
            updated_time=asset.created_at or datetime.datetime.now(),
            checksums=_archive_checksums(rgc, asset),
            access_methods=access_methods,
            contents=contents,
            description=f"Refgenie asset: {asset.registry_path}",
            aliases=[asset.digest],
        )

    # Archive only: a simple DRS object
    if "archive" in modes:
        access_methods = _build_access_methods(rgc, asset, base_uri, object_id, modes)

        return DRSModel(
            id=object_id,
            name=asset.registry_path,
            self_uri=f"drs://{request.url.netloc}/{object_id}",
            size=asset.size or 0,
            created_time=asset.created_at or datetime.datetime.now(),
            updated_time=asset.created_at or datetime.datetime.now(),
            checksums=_archive_checksums(rgc, asset),
            access_methods=access_methods,
            contents=None,
            description=f"Refgenie asset: {asset.registry_path}",
            aliases=[asset.digest],
        )

    # No serving mode: metadata only
    return DRSModel(
        id=object_id,
        name=asset.registry_path,
        self_uri=f"drs://{request.url.netloc}/{object_id}",
        size=asset.size or 0,
        created_time=asset.created_at or datetime.datetime.now(),
        updated_time=asset.created_at or datetime.datetime.now(),
        checksums=_archive_checksums(rgc, asset),
        access_methods=[],
        contents=None,
        description=f"Refgenie asset: {asset.registry_path}",
        aliases=[asset.digest],
    )


@router.get("/objects/{object_id}/access/{access_id}", response_model=AccessURL)
async def get_access_url(
    object_id: str,
    access_id: str,
    request: Request,
    rgc: Refgenie = Depends(get_refgenie),
):
    """Return a URL to fetch a DrsObject's bytes.

    ``archive``/``file`` give this server's ``/bytes`` URL;
    ``{type}:{remote_id}:{mode}`` gives the remote's direct URL.
    """
    asset_digest, file_path = parse_refgenie_drs_object_id(object_id)
    base_uri = f"{request.url.scheme}://{request.url.netloc}"

    if access_id in ("archive", "file"):
        if access_id == "file" and file_path is None:
            raise HTTPException(
                status_code=400,
                detail="file access requires file-level object ID ({digest}:{filename})",
            )
        require_published(rgc, asset_digest, access_id, file_path)
        return AccessURL(
            url=f"{base_uri}/ga4gh/drs/v1/objects/{object_id}/access/{access_id}/bytes",
            headers=None,
        )

    # Report a missing asset before any access_id parsing error.
    if not get_asset(rgc, asset_digest):
        raise HTTPException(status_code=404, detail=f"Asset not found: {asset_digest}")

    if ":" in access_id:
        parts = access_id.split(":")
        if len(parts) != 3:
            raise HTTPException(
                status_code=400,
                detail=f"Malformed access_id: {access_id}",
            )
        remote_type, remote_id_str, mode = parts
        try:
            remote_id = int(remote_id_str)
        except ValueError:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid remote_id in access_id: {access_id}",
            )

        asset, _ = require_published(rgc, asset_digest, mode, file_path)

        config = require_config(rgc)
        local_path = staged_local_path(config, asset, mode)
        if file_path is not None:
            local_path = local_path / file_path
        url = next(
            (
                u
                for remote, _, u in rgc.remote.pushed_urls(asset_digest, {mode: local_path})
                if remote.id == remote_id
            ),
            None,
        )
        if not url:
            raise HTTPException(
                status_code=404,
                detail=f"Remote URL not found for access_id: {access_id}",
            )
        return AccessURL(url=url, headers=None)

    raise HTTPException(status_code=404, detail=f"Unknown access method: {access_id}")


@router.get("/objects/{object_id}/access/{access_id}/bytes")
async def get_object_bytes(
    object_id: str,
    access_id: str,
    rgc: Refgenie = Depends(get_refgenie),
):
    """Return a DrsObject's bytes for the local access methods (``archive``, ``file``)."""
    asset_digest, file_path = parse_refgenie_drs_object_id(object_id)

    if access_id not in ("archive", "file"):
        raise HTTPException(
            status_code=400,
            detail=f"Only local access methods use /bytes: {access_id}",
        )

    if access_id == "archive":
        asset, _ = require_published(rgc, asset_digest, "archive")
        config = require_config(rgc)
        tarball_path = staged_local_path(config, asset, "archive")
        if not tarball_path.is_file():
            raise HTTPException(
                status_code=404,
                detail=f"Archive not found on disk for {asset_digest}",
            )
        return FileResponse(
            tarball_path,
            filename=tarball_path.name,
            media_type="application/gzip",
        )

    if access_id == "file":
        # Redirect to the v4 file-serving endpoint
        if file_path is None:
            raise HTTPException(
                status_code=400,
                detail="file access requires file-level object ID ({digest}:{filename})",
            )
        require_published(rgc, asset_digest, "file", file_path)
        return RedirectResponse(url=f"/v4/assets/{asset_digest}/files/{file_path}")
