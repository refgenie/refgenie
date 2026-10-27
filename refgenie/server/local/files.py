"""Read-only access to an asset's files on the user's disk (``/v1/files``).

Local mode only. The public server serves files through
``/v4/assets/{digest}/files/{path}``, which needs a *staged* (published) copy;
a local install has none, so the local dash reads the asset's own content
directory instead -- ``genome_folder / asset.path``, as recorded in the
database.

=======  ================================================  ====================
Method   Path                                              Result
=======  ================================================  ====================
GET      /v1/files/assets/{asset_digest}                   200 LocalAssetFiles
GET      /v1/files/assets/{asset_digest}/raw/{path}        200 bytes, attachment
GET      /v1/files/assets/{asset_digest}/preview/{path}    200 LocalFilePreview
=======  ================================================  ====================

Security
--------

* **Origin.** Every route depends on
  :func:`refgenie.server.local.security.require_local_origin`: the local SPA
  and the trusted dev origins only. The localhost bridge (the public site)
  gets 403 under every bridge mode -- it may learn *which* assets you have,
  never their bytes. Responses also carry ``Cross-Origin-Resource-Policy:
  same-origin`` and ``X-Content-Type-Options: nosniff`` so another page cannot
  embed them as a script or image either.
* **Paths.** A file is named by asset digest plus a path relative to the
  asset's directory. Absolute paths, ``..`` segments, backslashes and NUL are
  refused before the filesystem is touched (400 ``invalid_path``). The
  survivor is resolved with every symlink followed and must land inside the
  asset's directory -- or inside the same genome's content store
  (``genome_folder/data/<genome_digest>``), which is where colocation links
  (a child asset's link to its parent's file) point. Anything else, e.g. a
  link to ``~/.ssh``, is 404 ``file_not_found``, the same answer as a missing
  file, and is left out of the listing.
* **Preview.** At most :data:`PREVIEW_LIMIT` bytes, and only for text: a file
  that is compressed (by name or by its gzip/bgzip magic bytes) or whose
  first bytes are not UTF-8 text is download-only (409 ``not_previewable``).
"""

import os
import stat
from pathlib import Path, PurePosixPath

from fastapi import APIRouter, Depends, HTTPException, Response
from fastapi.responses import FileResponse
from sqlmodel import Session

from refgenie.core import Refgenie
from refgenie.server.dependencies import get_db_session, get_refgenie
from refgenie.server.errors import ErrorCode
from refgenie.server.local.schemas import LocalAssetFiles, LocalFileEntry, LocalFilePreview
from refgenie.server.local.security import require_local_origin
from refgenie.server.routers.helpers import get_asset_with_group

__all__ = ["PREVIEW_LIMIT", "MAX_LIST_ENTRIES", "router"]

#: Bytes shown by a preview. Enough for a .fai, chrom sizes or a JSON/YAML
#: descriptor in full, and the head of anything bigger.
PREVIEW_LIMIT = 64 * 1024

#: Bytes sniffed per file when the listing decides ``previewable``.
_SNIFF_BYTES = 4096

#: A listing stops here (``truncated: true``). Index directories can hold many
#: files; the UI is a list, not a file manager.
MAX_LIST_ENTRIES = 2000

#: Compressed or binary formats that are download-only whatever their bytes say.
_BINARY_SUFFIXES = frozenset(
    {
        ".gz",
        ".bgz",
        ".bz2",
        ".xz",
        ".zst",
        ".zip",
        ".tar",
        ".bam",
        ".cram",
        ".bai",
        ".crai",
        ".tbi",
        ".csi",
        ".bw",
        ".bigwig",
        ".bb",
        ".bigbed",
        ".2bit",
    }
)

#: The gzip magic number; bgzip (``.bgz``, ``.fa.gz``) shares it.
_GZIP_MAGIC = b"\x1f\x8b"

#: Applied to every byte-carrying response from this router.
_SAFE_HEADERS = {
    "Cross-Origin-Resource-Policy": "same-origin",
    "X-Content-Type-Options": "nosniff",
    "Cache-Control": "no-store",
}

router = APIRouter(
    prefix="/files",
    tags=["Local files"],
    dependencies=[Depends(require_local_origin)],
)


def _error(status: int, code: ErrorCode, message: str) -> HTTPException:
    return HTTPException(status_code=status, detail={"code": str(code), "message": message})


def _file_not_found(file_path: str) -> HTTPException:
    return _error(404, ErrorCode.FILE_NOT_FOUND, f"No file '{file_path}' in this asset.")


def validate_relative_path(file_path: str) -> PurePosixPath:
    """The request path as a safe relative path, or 400 ``invalid_path``.

    Purely lexical; the containment check after symlink resolution is
    :func:`resolve_inside`'s job. Both run on every request.
    """
    if (
        not file_path
        or "\x00" in file_path
        or "\\" in file_path
        or file_path.startswith("/")
        or (len(file_path) > 1 and file_path[1] == ":")  # C:foo
    ):
        raise _error(400, ErrorCode.INVALID_PATH, f"Invalid file path: {file_path!r}")
    path = PurePosixPath(file_path)
    if path.is_absolute() or any(part in ("..", ".", "") for part in file_path.split("/")):
        raise _error(400, ErrorCode.INVALID_PATH, f"Invalid file path: {file_path!r}")
    return path


def _is_within(path: Path, root: Path) -> bool:
    return path == root or path.is_relative_to(root)


def resolve_inside(candidate: Path, roots: tuple[Path, ...]) -> Path | None:
    """``candidate`` with every symlink resolved, if it lands inside a root.

    ``roots`` must already be resolved. Returns None for a missing path, a
    symlink loop, or anything that escapes every root.
    """
    try:
        resolved = candidate.resolve(strict=True)
    except (OSError, RuntimeError):
        return None
    return resolved if any(_is_within(resolved, root) for root in roots) else None


class _AssetDir:
    """An asset's content directory and the roots its files may resolve into."""

    def __init__(self, directory: Path, roots: tuple[Path, ...]):
        self.directory = directory
        self.roots = roots


def _asset_dir(session: Session, rgc: Refgenie, asset_digest: str) -> _AssetDir:
    """Look the asset up by digest and derive its directory from the DB record."""
    asset = get_asset_with_group(session, asset_digest)
    if asset is None:
        raise _error(404, ErrorCode.ASSET_NOT_FOUND, f"Asset {asset_digest} not found.")
    recorded = asset.path
    if not recorded or recorded.startswith("/") or ".." in PurePosixPath(recorded).parts:
        # A recorded path is always relative (data/<genome>/<group>/<digest>);
        # anything else is a corrupt row and is not followed.
        raise _error(404, ErrorCode.FILE_NOT_FOUND, f"Asset {asset_digest} has no local files.")
    genome_folder = Path(rgc.genome_folder)
    try:
        genome_root = genome_folder.resolve(strict=True)
        directory = (genome_folder / recorded).resolve(strict=True)
    except (OSError, RuntimeError):
        raise _error(
            404, ErrorCode.FILE_NOT_FOUND, f"Asset {asset_digest} has no files on disk."
        ) from None
    if not directory.is_dir() or not _is_within(directory, genome_root):
        raise _error(404, ErrorCode.FILE_NOT_FOUND, f"Asset {asset_digest} has no files on disk.")
    roots: tuple[Path, ...] = (directory,)
    genome_digest = asset.asset_group.genome_digest if asset.asset_group else None
    if genome_digest:
        try:
            content_store = (genome_folder / "data" / genome_digest).resolve(strict=True)
        except (OSError, RuntimeError):
            content_store = None
        if content_store is not None and _is_within(directory, content_store):
            roots = (directory, content_store)
    return _AssetDir(directory, roots)


def _resolve_file(asset_dir: _AssetDir, file_path: str) -> Path:
    relative = validate_relative_path(file_path)
    resolved = resolve_inside(asset_dir.directory / relative, asset_dir.roots)
    if resolved is None or not resolved.is_file():
        raise _file_not_found(file_path)
    return resolved


def _looks_binary_name(name: str) -> bool:
    return any(name.lower().endswith(suffix) for suffix in _BINARY_SUFFIXES)


def decode_text(head: bytes, truncated: bool) -> str | None:
    """``head`` as text, or None when it is binary.

    NUL bytes or invalid UTF-8 mean binary. When the head was cut short, a
    multi-byte character split at the cut is dropped rather than counted as
    invalid.
    """
    if head.startswith(_GZIP_MAGIC) or b"\x00" in head:
        return None
    try:
        return head.decode("utf-8")
    except UnicodeDecodeError as exc:
        if truncated and exc.start >= len(head) - 3:
            try:
                return head[: exc.start].decode("utf-8")
            except UnicodeDecodeError:
                return None
        return None


def _is_previewable(path: Path, name: str, size: int) -> bool:
    if _looks_binary_name(name):
        return False
    if size == 0:
        return True
    try:
        with path.open("rb") as handle:
            head = handle.read(_SNIFF_BYTES)
    except OSError:
        return False
    return decode_text(head, truncated=size > len(head)) is not None


@router.get("/assets/{asset_digest}", response_model=LocalAssetFiles)
def list_local_asset_files(
    asset_digest: str,
    response: Response,
    session: Session = Depends(get_db_session),
    rgc: Refgenie = Depends(get_refgenie),
):
    """List every file and directory in an asset's content directory.

    Recursive, sorted by path, directories included as their own entries.
    Symlinked directories are listed but not descended into; entries that
    resolve outside the asset's content are omitted.
    """
    asset_dir = _asset_dir(session, rgc, asset_digest)
    entries: list[LocalFileEntry] = []
    truncated = False
    for current, dirnames, filenames in os.walk(asset_dir.directory, followlinks=False):
        dirnames.sort()
        base = Path(current)
        for name in sorted(dirnames + filenames):
            if len(entries) >= MAX_LIST_ENTRIES:
                truncated = True
                break
            candidate = base / name
            resolved = resolve_inside(candidate, asset_dir.roots)
            if resolved is None:
                continue
            relative = candidate.relative_to(asset_dir.directory).as_posix()
            try:
                info = resolved.stat()
            except OSError:
                continue
            if stat.S_ISDIR(info.st_mode):
                entries.append(LocalFileEntry(path=relative, type="dir"))
            elif stat.S_ISREG(info.st_mode):
                entries.append(
                    LocalFileEntry(
                        path=relative,
                        type="file",
                        size=info.st_size,
                        previewable=_is_previewable(resolved, name, info.st_size),
                    )
                )
        if truncated:
            break
    entries.sort(key=lambda entry: entry.path)
    response.headers.update(_SAFE_HEADERS)
    return LocalAssetFiles(asset_digest=asset_digest, entries=entries, truncated=truncated)


@router.get("/assets/{asset_digest}/raw/{file_path:path}", response_class=FileResponse)
def download_local_asset_file(
    asset_digest: str,
    file_path: str,
    session: Session = Depends(get_db_session),
    rgc: Refgenie = Depends(get_refgenie),
):
    """Download one file of an asset (``Content-Disposition: attachment``)."""
    path = _resolve_file(_asset_dir(session, rgc, asset_digest), file_path)
    return FileResponse(
        path,
        filename=PurePosixPath(file_path).name,
        media_type="application/octet-stream",
        content_disposition_type="attachment",
        headers=_SAFE_HEADERS,
    )


@router.get(
    "/assets/{asset_digest}/preview/{file_path:path}",
    response_model=LocalFilePreview,
)
def preview_local_asset_file(
    asset_digest: str,
    file_path: str,
    response: Response,
    session: Session = Depends(get_db_session),
    rgc: Refgenie = Depends(get_refgenie),
):
    """The first :data:`PREVIEW_LIMIT` bytes of a text file, as JSON.

    409 ``not_previewable`` for compressed or binary files.
    """
    path = _resolve_file(_asset_dir(session, rgc, asset_digest), file_path)
    not_previewable = _error(
        409,
        ErrorCode.NOT_PREVIEWABLE,
        f"'{file_path}' is binary or compressed; download it instead.",
    )
    if _looks_binary_name(PurePosixPath(file_path).name) or _looks_binary_name(path.name):
        raise not_previewable
    size = path.stat().st_size
    with path.open("rb") as handle:
        head = handle.read(PREVIEW_LIMIT)
    truncated = size > len(head)
    content = decode_text(head, truncated)
    if content is None:
        raise not_previewable
    response.headers.update(_SAFE_HEADERS)
    return LocalFilePreview(
        path=file_path, size=size, content=content, truncated=truncated, limit=PREVIEW_LIMIT
    )
