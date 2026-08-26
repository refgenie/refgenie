"""One error vocabulary and one classifier for the whole web layer.

Every failure in local mode reports the same `code` whether it surfaces as an
HTTP error body or as a `JobError` on a background job record. That is the
point of this module: the classifier lives here, is consumed by the HTTP
exception handlers *and* by `refgenie/server/jobs/runners.py`, and is defined
in exactly one place so the two can never drift.

Imports nothing from `refgenie.server`, so both consumers can import it without
a cycle.

Ordering matters. `AssetExistsError` and `NoArchiveError` both subclass
`PullFailedError`, so `classify_exception` scans an **ordered** list with
`isinstance`, most specific first. A type-keyed dict or a mis-ordered
`except` chain silently collapses the three into one, which is the exact
differentiation the pull UI exists to show. `tests/test_db.py` pins the
ordering.
"""

import logging
from enum import StrEnum
from typing import TYPE_CHECKING

from pydantic import BaseModel

from refgenie.exceptions import (
    AssetClassExistsError,
    AssetExistsError,
    MissingAliasError,
    MissingAssetClassError,
    MissingAssetError,
    MissingAssetGroupError,
    MissingBuildInputError,
    MissingGenomeError,
    MissingRecipeError,
    NoArchiveError,
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
    """The single machine-readable error vocabulary.

    Referenced by both `ErrorBody.code` (HTTP) and `JobError.code` (jobs). The
    UI branches on these, never on a message string.
    """

    # --- pull -------------------------------------------------------------
    ASSET_EXISTS = "asset_exists"
    NO_ARCHIVE = "no_archive"
    PULL_FAILED = "pull_failed"
    PULL_SKIPPED = "pull_skipped"
    SERVER_CANNOT_SERVE = "server_cannot_serve"
    REMOTE_DIGEST_MISMATCH = "remote_digest_mismatch"
    NO_SUBSCRIPTIONS = "no_subscriptions"

    # --- build ------------------------------------------------------------
    BUILD_FAILED = "build_failed"
    MISSING_BUILD_INPUT = "missing_build_input"

    # --- jobs -------------------------------------------------------------
    CANCELLED = "cancelled"

    # --- lookups ----------------------------------------------------------
    #: Generic 404 for HTTPExceptions raised with a plain-string detail (a job
    #: id, a page). The resource-specific *_not_found codes below are what the
    #: actions endpoints emit; this one only backfills the envelope for routes
    #: that predate it.
    NOT_FOUND = "not_found"
    GENOME_NOT_FOUND = "genome_not_found"
    ALIAS_NOT_FOUND = "alias_not_found"
    ASSET_NOT_FOUND = "asset_not_found"
    ASSET_GROUP_NOT_FOUND = "asset_group_not_found"
    ASSET_CLASS_NOT_FOUND = "asset_class_not_found"
    RECIPE_NOT_FOUND = "recipe_not_found"

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


class WebError(Exception):
    """A failure the web layer synthesizes rather than catches.

    Some outcomes are not exceptions in the library. `Refgenie.pull` returns
    `None` when nothing is subscribed and `build_asset` returns `None` when the
    pipeline failed -- correct for a CLI that prints and exits, useless for a
    job record that must say what went wrong. Raising one of these carries the
    right code through the same classifier as everything else, instead of
    growing a second mapping for "things that were not exceptions".
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


#: Codes whose user-facing message is written here rather than taken from the
#: exception. `AssetExistsError`'s own text names an internal registry path;
#: "use force" is the actionable half and the UI shows it verbatim.
_MESSAGE_OVERRIDES: dict[ErrorCode, str] = {
    ErrorCode.ASSET_EXISTS: "Asset already exists. Use force to overwrite.",
}

#: (exception type, HTTP status, code), scanned in order with isinstance.
#:
#: MOST SPECIFIC FIRST. AssetExistsError and NoArchiveError subclass
#: PullFailedError; MissingAliasError and the rest subclass RefgenieError.
#: Reordering this list changes behavior -- it is not a stylistic list.
_CLASSIFICATION: tuple[tuple[type[BaseException], int, ErrorCode], ...] = (
    # PullFailedError subclasses, before PullFailedError itself
    (AssetExistsError, 409, ErrorCode.ASSET_EXISTS),
    (NoArchiveError, 404, ErrorCode.NO_ARCHIVE),
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
    # Catch-alls, last
    (RefgenieError, 500, ErrorCode.REFGENIE_ERROR),
    (ValueError, 409, ErrorCode.CONFLICT),
)


def classify_exception(exc: BaseException) -> tuple[int, ErrorCode, str]:
    """Map an exception to (HTTP status, error code, message).

    The one place a refgenie exception becomes a machine-readable code. Used by
    the local-mode HTTP handlers and by the job runners, so `asset_exists`
    means the same thing on both transports.

    Args:
        exc: The exception raised by the facade.

    Returns:
        A (status, code, message) triple. Unrecognized exceptions get
        (500, `internal_error`, str(exc)).
    """
    # A synthesized error already knows its code; it is not looked up.
    if isinstance(exc, WebError):
        return exc.status, exc.code, str(exc)
    for exc_type, status, code in _CLASSIFICATION:
        if isinstance(exc, exc_type):
            return status, code, _MESSAGE_OVERRIDES.get(code, str(exc))
    return 500, ErrorCode.INTERNAL_ERROR, str(exc)


# ---------------------------------------------------------------------------
# The HTTP envelope
# ---------------------------------------------------------------------------
#
# Everything below is HTTP-side sugar over the classifier. It stays in this
# module so there is exactly one file where a code is minted, but the fastapi
# imports are deferred into the functions: `jobs/schemas.py` imports this module
# for `ErrorCode` alone and must stay import-cheap.


def error_response(
    status: int,
    code: "ErrorCode | str",
    message: str,
    detail: str | None = None,
    field: str | None = None,
):
    """A JSONResponse carrying the standard envelope.

    The body is built as a plain dict rather than through `ErrorBody` on
    purpose: a few router-minted codes (`not_cancellable`, `already_done`)
    are legitimate strings outside the enum, and the envelope must be able to
    carry them.
    """
    from fastapi.responses import JSONResponse

    error: dict = {"code": str(code), "message": message}
    if detail is not None:
        error["detail"] = detail
    if field is not None:
        error["field"] = field
    return JSONResponse(status_code=status, content={"ok": False, "error": error})


def action_error(exc: BaseException) -> "HTTPException":
    """The HTTPException for a classified failure, envelope detail attached.

    Handlers raise this from their narrow `except` clauses::

        try:
            rgc.alias.remove(alias_name)
        except MissingAliasError as exc:
            raise action_error(exc) from exc
    """
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
    """Make every non-2xx body from this app use the one envelope.

    Local mode only. Three handlers:

    * ``HTTPException`` -- a dict detail that already looks like an error body
      (has ``code`` and ``message``) passes through unchanged; anything else is
      wrapped with a code derived from the status.
    * ``RequestValidationError`` -- 422, code ``validation_error``, with the
      first offending field named.
    * ``Exception`` -- 500, code ``internal_error``. The traceback goes to the
      log, never to the browser (the old dash handler returned the *handler's*
      stack to the client; do not resurrect that).
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
