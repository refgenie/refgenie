import datetime
import logging

from fastapi import APIRouter, Depends, HTTPException, Request
from fastapi.responses import FileResponse
from refgenie.db.tables import (
    Asset,
    AssetGroup,
    Remote,
    RemoteAssetLink,
)
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select
from starlette.responses import RedirectResponse

from refgenie.server.dependencies import get_db_session
from refgenie.server.remote_assets import get_remote_url
from refgenie.server.routers.helpers import (
    get_asset_with_group,
    get_config,
    get_staged_asset,
    require_config,
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
    """
    Parse a Refgenie DRS object ID into (asset_digest, file_path).

    Two forms:
      - "{64-hex-digest}" -> (digest, None)        -- asset-level
      - "{64-hex-digest}:{filename}" -> (digest, filename) -- file-level

    Colon is URL-safe in path segments so this is valid in DRS object IDs.
    """
    if not object_id or len(object_id) < 8:
        raise HTTPException(status_code=400, detail=f"Object ID {object_id} is malformed")

    if ":" in object_id:
        digest, file_path = object_id.split(":", 1)
        return digest, file_path
    return object_id, None


def _get_serving_modes(asset: Asset) -> set[str]:
    """Resolve effective serving modes for an asset."""
    return set(asset.serving_modes)


def _archive_checksums(session: Session, asset: Asset) -> list[Checksum]:
    """
    Build DRS ``checksums`` describing the bytes a client actually downloads.

    A DRS checksum must match the bytes served through an access method, NOT
    the asset's content-identity digest. For archive serving, the downloaded
    bytes are the tarball, whose sha-256 is ``StagedAsset.tarball_digest``.
    (``asset.digest`` is the content-identity digest -- correct for the DRS
    object id, aliases, self_uri, and S3 key, but it is NOT the checksum of
    the tarball bytes.)

    Returns an empty list when no archive tarball digest is available (e.g.
    file-only or metadata-only assets) rather than advertising a wrong
    checksum. This selects the ``mode == "archive"`` staged asset explicitly.
    """
    staged_archive = get_staged_asset(session, asset.digest, "archive")
    if staged_archive and staged_archive.tarball_digest:
        return [Checksum(checksum=staged_archive.tarball_digest, type="sha-256")]
    return []


def _build_access_methods(
    session: Session,
    asset: Asset,
    base_uri: str,
    object_id: str,
    modes: set[str],
) -> list[AccessMethod]:
    """
    Build DRS access methods for all available locations.

    For each serving mode, adds:
    - A local access method (if genome_stage_folder is configured and staged asset exists)

    access_id scheme:
    - "archive"                              -- local archive (this server)
    - "file"                                 -- local files (this server)
    - "{remote_type}:{remote_id}:archive"    -- remote archive (e.g., "s3:1:archive")
    - "{remote_type}:{remote_id}:file"       -- remote files (e.g., "https:2:file")
    """
    methods = []

    # 1. Local access methods (only if server has genome_stage_folder on disk)
    config = get_config(session)
    has_local = config and config.genome_stage_folder

    if has_local:
        if "archive" in modes:
            staged_archive = get_staged_asset(session, asset.digest, "archive")
            if staged_archive:
                methods.append(
                    AccessMethod(
                        type="https",
                        access_url=AccessURL(
                            url=f"{base_uri}/ga4gh/drs/objects/{object_id}/access/archive",
                            headers=None,
                        ),
                        access_id="archive",
                    )
                )

        if "file" in modes:
            staged_file = get_staged_asset(session, asset.digest, "file")
            if staged_file:
                methods.append(
                    AccessMethod(
                        type="https",
                        access_url=AccessURL(
                            url=f"{base_uri}/ga4gh/drs/objects/{object_id}/access/file",
                            headers=None,
                        ),
                        access_id="file",
                    )
                )

    # 2. Remote access methods (from pushed RemoteAssetLinks)
    remote_links = session.exec(
        select(RemoteAssetLink, Remote)
        .join(Remote, RemoteAssetLink.remote_id == Remote.id)
        .where(RemoteAssetLink.asset_digest == asset.digest)
        .where(RemoteAssetLink.pushed == True)  # noqa: E712
        .where(RemoteAssetLink.mode.in_(list(modes)))
    ).all()

    # Derive a local path to compute the relative path for remote URL
    asset_with_group = get_asset_with_group(session, asset.digest) if remote_links else None

    for link, remote in remote_links:
        access_id_str = f"{remote.type.value}:{remote.id}:{link.mode}"
        if asset_with_group and has_local:
            local_path = staged_local_path(config, asset_with_group, link.mode)
            remote_url = get_remote_url(local_path, session, asset.digest, link.mode)
            if remote_url:
                methods.append(
                    AccessMethod(
                        type=remote.type.value,
                        access_url=AccessURL(url=remote_url, headers=None),
                        access_id=access_id_str,
                    )
                )

    return methods


@router.get("/service-info", response_model=ServiceInfo)
async def get_service_info(request: Request):
    """
    Returns information about the DRS service
    """
    return ServiceInfo(
        id="org.refgenie.server",
        name="Refgenie Server DRS",
        type=ServiceType(group="org.ga4gh", artifact="drs", version=DRS_SPEC_VERSION),
        description=(
            "Refgenie Server implementing GA4GH DRS specification for genomic reference assets"
        ),
        organization=Organization(name="Refgenie", url="https://refgenie.databio.org"),
        contactUrl="https://github.com/refgenie/refgenieserver/issues",
        documentationUrl="https://refgenie.databio.org",
        createdAt=datetime.datetime.now(),
        updatedAt=datetime.datetime.now(),
        environment="production",
        version=DRS_SPEC_VERSION,
    )


@router.get("/objects/{object_id}", response_model=DRSModel)
async def get_drs_object_metadata(
    object_id: str,
    request: Request,
    session: Session = Depends(get_db_session),
):
    """
    Returns metadata about a DrsObject.

    Serving-mode-aware: returns structurally different responses depending on
    how the asset is served (file, archive, both, or none).

    For file-level requests ({digest}:{filename}), returns a simple DRS object
    with file-mode access methods only.
    """
    asset_digest, file_path = parse_refgenie_drs_object_id(object_id)
    base_uri = f"{request.url.scheme}://{request.url.netloc}"

    with session:
        # Eager-load asset_group -> asset_class so serving_modes resolves without lazy loading
        asset = (
            session.exec(
                select(Asset)
                .options(selectinload(Asset.asset_group).selectinload(AssetGroup.asset_class))
                .where(Asset.digest == asset_digest)
            )
            .unique()
            .one_or_none()
        )

        if not asset:
            raise HTTPException(
                status_code=404, detail=f"Asset with digest {asset_digest} not found"
            )

        modes = _get_serving_modes(asset)

        # --- File-level request ({digest}:{filename}) ---
        if file_path is not None:
            if "file" not in modes:
                raise HTTPException(
                    status_code=404,
                    detail=f"Asset {asset_digest} does not support file-level access",
                )

            # Verify file exists in directory_contents
            staged_file = get_staged_asset(session, asset.digest, "file")

            if not staged_file or file_path not in (staged_file.directory_contents or []):
                raise HTTPException(
                    status_code=404,
                    detail=f"File '{file_path}' not found in asset {asset_digest}",
                )

            access_methods = _build_access_methods(
                session, asset, base_uri, object_id, modes={"file"}
            )

            return DRSModel(
                id=object_id,
                name=file_path,
                self_uri=f"drs://{request.url.netloc}/{object_id}",
                size=0,  # Individual file sizes not tracked yet; 0 is spec-legal
                created_time=asset.created_at or datetime.datetime.now(),
                updated_time=asset.created_at or datetime.datetime.now(),
                # File-level object serves a single file via file mode. There is
                # no tracked per-file byte digest, so omit rather than advertise a
                # wrong sha-256 (asset.digest is the identity digest, not the file's).
                checksums=[],
                access_methods=access_methods,
                contents=None,  # Simple object, not a bundle
                description=f"File: {file_path}",
                aliases=[],
            )

        # --- Asset-level request ---

        # Case: file mode (with or without archive) -> DRS bundle
        if "file" in modes:
            staged_file = get_staged_asset(session, asset.digest, "file")

            contents = [
                ContentsObject(
                    name=filename,
                    id=f"{asset.digest}:{filename}",
                    drs_uri=[f"drs://{request.url.netloc}/{asset.digest}:{filename}"],
                    contents=None,
                )
                for filename in (staged_file.directory_contents if staged_file else [])
            ]

            access_methods = _build_access_methods(session, asset, base_uri, object_id, modes)

            return DRSModel(
                id=object_id,
                name=asset.registry_path,
                self_uri=f"drs://{request.url.netloc}/{object_id}",
                size=asset.size or 0,
                created_time=asset.created_at or datetime.datetime.now(),
                updated_time=asset.created_at or datetime.datetime.now(),
                # Advertise the archive tarball digest when an archive mode is
                # available (combo assets); pure file-mode bundles have no single
                # downloadable-byte digest, so this is empty for them.
                checksums=_archive_checksums(session, asset),
                access_methods=access_methods,
                contents=contents,
                description=f"Refgenie asset: {asset.registry_path}",
                aliases=[asset.digest],
            )

        # Case: archive only (no file) -> simple DRS object
        if "archive" in modes:
            access_methods = _build_access_methods(session, asset, base_uri, object_id, modes)

            return DRSModel(
                id=object_id,
                name=asset.registry_path,
                self_uri=f"drs://{request.url.netloc}/{object_id}",
                size=asset.size or 0,
                created_time=asset.created_at or datetime.datetime.now(),
                updated_time=asset.created_at or datetime.datetime.now(),
                # Archive-mode object: the client downloads the tarball, so the
                # checksum is the tarball byte digest, not the identity digest.
                checksums=_archive_checksums(session, asset),
                access_methods=access_methods,
                contents=None,
                description=f"Refgenie asset: {asset.registry_path}",
                aliases=[asset.digest],
            )

        # Case: none mode -> metadata only
        return DRSModel(
            id=object_id,
            name=asset.registry_path,
            self_uri=f"drs://{request.url.netloc}/{object_id}",
            size=asset.size or 0,
            created_time=asset.created_at or datetime.datetime.now(),
            updated_time=asset.created_at or datetime.datetime.now(),
            # None mode serves no bytes and has no tarball, so no checksum.
            checksums=_archive_checksums(session, asset),
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
    session: Session = Depends(get_db_session),
):
    """
    Returns a URL that can be used to fetch the bytes of a DrsObject.

    Routes by access_id:
    - "archive" or "file": local access, returns URL to /bytes endpoint
    - "{type}:{remote_id}:{mode}": remote access, returns direct remote URL
    """
    asset_digest, file_path = parse_refgenie_drs_object_id(object_id)
    base_uri = f"{request.url.scheme}://{request.url.netloc}"

    with session:
        asset = get_asset_with_group(session, asset_digest)
        if not asset:
            raise HTTPException(status_code=404, detail=f"Asset not found: {asset_digest}")

        # Parse access_id: "archive", "file", or "{type}:{remote_id}:{mode}"
        if access_id in ("archive", "file"):
            # Local access -- serve from this server's /bytes endpoint
            return AccessURL(
                url=f"{base_uri}/ga4gh/drs/objects/{object_id}/access/{access_id}/bytes",
                headers=None,
            )

        if ":" in access_id:
            # Remote access -- parse and resolve
            parts = access_id.split(":")
            if len(parts) != 3:
                raise HTTPException(
                    status_code=400,
                    detail=f"Malformed access_id: {access_id}",
                )
            remote_type, remote_id_str, mode = parts
            try:
                int(remote_id_str)
            except ValueError:
                raise HTTPException(
                    status_code=400,
                    detail=f"Invalid remote_id in access_id: {access_id}",
                )

            # Resolve remote URL via get_remote_url
            config = require_config(session)
            local_path = staged_local_path(config, asset, mode)

            url = get_remote_url(local_path, session, asset_digest, mode)
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
    session: Session = Depends(get_db_session),
):
    """
    Returns the bytes of a DrsObject for local access methods.

    Only local access_ids ("archive", "file") reach this endpoint.
    Remote access methods resolve to direct URLs in get_access_url().
    """
    asset_digest, file_path = parse_refgenie_drs_object_id(object_id)

    # Only local access_ids reach /bytes
    if access_id not in ("archive", "file"):
        raise HTTPException(
            status_code=400,
            detail=f"Only local access methods use /bytes: {access_id}",
        )

    with session:
        asset = get_asset_with_group(session, asset_digest)
        if not asset:
            raise HTTPException(status_code=404, detail=f"Asset not found: {asset_digest}")

        config = require_config(session)

        if access_id == "archive":
            # Serve tarball from local disk (content-addressed)
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
            # Serve individual file -- redirect to the v4 file-serving endpoint
            if file_path is None:
                raise HTTPException(
                    status_code=400,
                    detail="file access requires file-level object ID ({digest}:{filename})",
                )
            return RedirectResponse(url=f"/v4/assets/{asset_digest}/files/{file_path}")
