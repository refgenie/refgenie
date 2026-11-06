"""The web layer's error codes and exception classifier.

`classify_exception` is used by the HTTP handlers and by the job manager, so a
failure reports the same `code` in an HTTP body and on a `JobError`. This
module imports nothing from `refgenie.server`, to avoid import cycles.
"""

import logging
from enum import StrEnum
from typing import TYPE_CHECKING

from pydantic import BaseModel

from refgenie.exceptions import (
    AssetClassExistsError,
    AssetExistsError,
    BuildFailedError,
    ContainerEngineMissingError,
    MissingAliasError,
    MissingAssetClassError,
    MissingAssetError,
    MissingAssetGroupError,
    MissingBuildInputError,
    MissingGenomeError,
    MissingRecipeError,
    NoArchiveError,
    NoSubscriptionsError,
    PullFailedError,
    PullSkipped,
    RecipeExistsError,
    RefgenieError,
    RemoteDigestMismatchError,
    ServerCannotServe,
)

if TYPE_CHECKING:
    from fastapi import FastAPI, HTTPException

__all__ = [
    "ErrorBody",
    "ErrorCode",
    "ErrorResponse",
    "WebError",
    "action_error",
    "classify_exception",
    "error_response",
    "install_error_handlers",
]

_log = logging.getLogger(__name__)


class ErrorCode(StrEnum):
    """Machine-readable error codes, used by `ErrorBody` and `JobError`."""

    # --- pull -------------------------------------------------------------
    ASSET_EXISTS = "asset_exists"
    NO_ARCHIVE = "no_archive"
    PULL_FAILED = "pull_failed"
    PULL_SKIPPED = "pull_skipped"
    SERVER_CANNOT_SERVE = "server_cannot_serve"
    REMOTE_DIGEST_MISMATCH = "remote_digest_mismatch"
    NO_SUBSCRIPTIONS = "no_subscriptions"
    #: ``server_url`` names a server this instance does not subscribe to.
    SERVER_NOT_SUBSCRIBED = "server_not_subscribed"
    #: A remote listing failed. The cause is logged, not returned.
    REMOTE_UNAVAILABLE = "remote_unavailable"

    # --- build ------------------------------------------------------------
    BUILD_FAILED = "build_failed"
    MISSING_BUILD_INPUT = "missing_build_input"
    #: Preflight only: the recipe names its output after a tool version that
    #: could not be probed, so the caller must supply ``asset``.
    ASSET_NAME_REQUIRED = "asset_name_required"
    #: The recipe needs a container and this host has neither Docker nor Apptainer.
    CONTAINER_ENGINE_MISSING = "container_engine_missing"

    # --- jobs -------------------------------------------------------------
    CANCELLED = "cancelled"

    # --- lookups ----------------------------------------------------------
    #: Generic 404 for an HTTPException with a plain-string detail.
    NOT_FOUND = "not_found"
    GENOME_NOT_FOUND = "genome_not_found"
    #: The genome is missing but aliases still name its digest; remove the alias.
    STALE_ALIAS = "stale_alias"
    ALIAS_NOT_FOUND = "alias_not_found"
    ASSET_NOT_FOUND = "asset_not_found"
    ASSET_GROUP_NOT_FOUND = "asset_group_not_found"
    ASSET_CLASS_NOT_FOUND = "asset_class_not_found"
    RECIPE_NOT_FOUND = "recipe_not_found"
    #: ``/v1/files``: the path is not a file inside the asset (missing, a
    #: directory, or a symlink that leads out of the asset's content).
    FILE_NOT_FOUND = "file_not_found"

    # --- generic ----------------------------------------------------------
    ALREADY_EXISTS = "already_exists"
    CONFLICT = "conflict"
    REFGENIE_ERROR = "refgenie_error"
    INTERNAL_ERROR = "internal_error"

    # --- router-only ------------------------------------------------------
    MISSING_ACTION_HEADER = "missing_action_header"
    FORBIDDEN_HOST = "forbidden_host"
    FORBIDDEN_ORIGIN = "forbidden_origin"
    VALIDATION_ERROR = "validation_error"
    #: ``/v1/files``: an absolute path, a ``..`` segment, or another path the
    #: route refuses before touching the filesystem.
    INVALID_PATH = "invalid_path"
    #: ``/v1/files``: the file is binary or compressed; download it instead.
    NOT_PREVIEWABLE = "not_previewable"


class WebError(Exception):
    """An error with an explicit code, for outcomes the library does not raise.

    For example, `build.run` returns `None` on failure; the build runner raises
    this instead.
    """

    def __init__(self, message: str, code: ErrorCode, status: int = 500):
        super().__init__(message)
        self.code = code
        self.status = status


class ErrorBody(BaseModel):
    """The `error` member of every non-2xx body from local mode."""

    code: ErrorCode
    message: str
    detail: str | None = None
    field: str | None = None


class ErrorResponse(BaseModel):
    """The envelope: `{"ok": false, "error": {...}}`."""

    ok: bool = False
    error: ErrorBody


#: Codes whose user-facing message replaces the exception's text.
_MESSAGE_OVERRIDES: dict[ErrorCode, str] = {
    ErrorCode.ASSET_EXISTS: "Asset already exists. Use force to overwrite.",
}

#: (exception type, HTTP status, code), scanned in order with isinstance, so
#: subclasses must come before their parents.
_CLASSIFICATION: tuple[tuple[type[BaseException], int, ErrorCode], ...] = (
    # PullFailedError subclasses, before PullFailedError itself
    (AssetExistsError, 409, ErrorCode.ASSET_EXISTS),
    (NoArchiveError, 404, ErrorCode.NO_ARCHIVE),
    (NoSubscriptionsError, 409, ErrorCode.NO_SUBSCRIPTIONS),
    (PullFailedError, 502, ErrorCode.PULL_FAILED),
    (PullSkipped, 409, ErrorCode.PULL_SKIPPED),
    (ServerCannotServe, 502, ErrorCode.SERVER_CANNOT_SERVE),
    (RemoteDigestMismatchError, 409, ErrorCode.REMOTE_DIGEST_MISMATCH),
    # Lookups
    (MissingGenomeError, 404, ErrorCode.GENOME_NOT_FOUND),
    (MissingAliasError, 404, ErrorCode.ALIAS_NOT_FOUND),
    (MissingAssetError, 404, ErrorCode.ASSET_NOT_FOUND),
    (MissingAssetGroupError, 404, ErrorCode.ASSET_GROUP_NOT_FOUND),
    (MissingAssetClassError, 404, ErrorCode.ASSET_CLASS_NOT_FOUND),
    (MissingRecipeError, 404, ErrorCode.RECIPE_NOT_FOUND),
    # Conflicts
    (RecipeExistsError, 409, ErrorCode.ALREADY_EXISTS),
    (AssetClassExistsError, 409, ErrorCode.ALREADY_EXISTS),
    (MissingBuildInputError, 400, ErrorCode.MISSING_BUILD_INPUT),
    (ContainerEngineMissingError, 400, ErrorCode.CONTAINER_ENGINE_MISSING),
    (BuildFailedError, 500, ErrorCode.BUILD_FAILED),
    # Catch-alls, last
    (RefgenieError, 500, ErrorCode.REFGENIE_ERROR),
    (ValueError, 409, ErrorCode.CONFLICT),
)


def classify_exception(exc: BaseException) -> tuple[int, ErrorCode, str]:
    """Map an exception to (HTTP status, error code, message).

    Args:
        exc: The exception raised by the root.

    Returns:
        A (status, code, message) triple. Unrecognized exceptions get
        (500, `internal_error`, str(exc)).
    """
    if isinstance(exc, WebError):
        return exc.status, exc.code, str(exc)
    for exc_type, status, code in _CLASSIFICATION:
        if isinstance(exc, exc_type):
            return status, code, _MESSAGE_OVERRIDES.get(code, str(exc))
    return 500, ErrorCode.INTERNAL_ERROR, str(exc)


# ---------------------------------------------------------------------------
# The HTTP envelope (fastapi imported lazily to keep this module cheap)
# ---------------------------------------------------------------------------


def error_response(
    status: int,
    code: "ErrorCode | str",
    message: str,
    detail: str | None = None,
    field: str | None = None,
):
    """A JSONResponse carrying the standard envelope.

    Built as a plain dict so codes outside `ErrorCode` (e.g. `not_cancellable`)
    can be carried.
    """
    from fastapi.responses import JSONResponse

    error: dict = {"code": str(code), "message": message}
    if detail is not None:
        error["detail"] = detail
    if field is not None:
        error["field"] = field
    return JSONResponse(status_code=status, content={"ok": False, "error": error})


def action_error(exc: BaseException) -> "HTTPException":
    """The HTTPException for a classified failure, envelope detail attached."""
    from fastapi import HTTPException

    status, code, message = classify_exception(exc)
    return HTTPException(status_code=status, detail={"code": str(code), "message": message})


def _code_for_status(status: int) -> str:
    """The fallback code for an HTTPException raised with a plain-string detail."""
    if status == 404:
        return str(ErrorCode.NOT_FOUND)
    if status == 409:
        return str(ErrorCode.CONFLICT)
    if status in (400, 422):
        return str(ErrorCode.VALIDATION_ERROR)
    if status == 421:
        return str(ErrorCode.FORBIDDEN_HOST)
    if status >= 500:
        return str(ErrorCode.INTERNAL_ERROR)
    return str(ErrorCode.REFGENIE_ERROR)


def install_error_handlers(app: "FastAPI") -> None:
    """Make every non-2xx body from this app (local mode) use the one envelope.

    HTTPException details with ``code`` and ``message`` pass through; others get
    a code from the status. Validation errors name the first bad field.
    Unhandled exceptions return a bare 500; the traceback goes only to the log.
    """
    from fastapi.exceptions import RequestValidationError
    from starlette.exceptions import HTTPException as StarletteHTTPException

    @app.exception_handler(StarletteHTTPException)
    async def _http_exception_handler(request, exc: StarletteHTTPException):
        from fastapi.responses import JSONResponse

        detail = exc.detail
        if isinstance(detail, dict) and "code" in detail and "message" in detail:
            error = detail
        else:
            error = {"code": _code_for_status(exc.status_code), "message": str(detail)}
        return JSONResponse(
            status_code=exc.status_code,
            content={"ok": False, "error": error},
            headers=getattr(exc, "headers", None),
        )

    @app.exception_handler(RequestValidationError)
    async def _validation_handler(request, exc: RequestValidationError):
        errors = exc.errors()
        first = errors[0] if errors else {}
        field = ".".join(str(part) for part in first.get("loc", ()) if part != "body") or None
        return error_response(
            422,
            ErrorCode.VALIDATION_ERROR,
            first.get("msg", "Invalid request body."),
            detail=str(errors),
            field=field,
        )

    @app.exception_handler(Exception)
    async def _unhandled_handler(request, exc: Exception):
        _log.exception(f"Unhandled error on {request.method} {request.url.path}")
        return error_response(500, ErrorCode.INTERNAL_ERROR, "Internal error")
