"""Request models and response envelopes for the local actions API.

Separate from the CLI command models, which key their schemas by CLI aliases.
Request models forbid unknown fields, so a misspelled field is a 422. The
pull, build and genome-init bodies are the job parameter models in
``refgenie.server.jobs.schemas``.
"""

from typing import Any, Literal
from urllib.parse import urlsplit

from pydantic import BaseModel, ConfigDict, Field, model_validator

__all__ = [
    "ActionRequest",
    "ActionResult",
    "AliasSetRequest",
    "DataChannelAddRequest",
    "LocalAssetFiles",
    "LocalFileEntry",
    "LocalFilePreview",
    "PreflightError",
    "PreflightResult",
    "SetDefaultAssetRequest",
    "SubscribeRequest",
    "UnsubscribeRequest",
]


class ActionRequest(BaseModel):
    """Base for every request model: unknown fields are an error."""

    model_config = ConfigDict(extra="forbid")


class AliasSetRequest(ActionRequest):
    """Body of ``POST /v1/actions/aliases``. The digest is required, to avoid remote lookups."""

    alias: str
    genome_digest: str


class SubscribeRequest(ActionRequest):
    """Body of ``POST /v1/actions/subscriptions``."""

    server_urls: list[str] = Field(min_length=1)
    reset: bool = False


class UnsubscribeRequest(ActionRequest):
    """Body of ``DELETE /v1/actions/subscriptions``."""

    server_urls: list[str] = Field(min_length=1)


class DataChannelAddRequest(ActionRequest):
    """Body of ``POST /v1/actions/data_channels``.

    Only ``http(s)`` indexes, since a local path would allow reading any local
    file. The channel type comes from the URL scheme.
    """

    name: str = Field(min_length=1)
    index_address: str
    description: str | None = None
    #: Register the channel's asset classes and recipes right away.
    sync: bool = True

    @model_validator(mode="after")
    def _http_index_only(self) -> "DataChannelAddRequest":
        scheme = urlsplit(self.index_address).scheme
        if scheme not in ("http", "https"):
            raise ValueError(
                "index_address must be an http:// or https:// URL to an index.yaml; "
                "local paths and ftp channels can only be added from the CLI."
            )
        return self

    @property
    def channel_type(self) -> str:
        return urlsplit(self.index_address).scheme


class SetDefaultAssetRequest(ActionRequest):
    """Body of ``POST /v1/actions/assets/default``."""

    genome_digest: str
    asset_group: str
    asset: str


# ---------------------------------------------------------------------------
# Response envelopes (a 202 returns the jobs package's ``JobRef``)
# ---------------------------------------------------------------------------


class ActionResult(BaseModel):
    """The 200 body of every synchronous action."""

    ok: Literal[True] = True
    message: str
    data: dict[str, Any] | None = None


class PreflightError(BaseModel):
    """One field-scoped problem found by a build preflight."""

    field: str
    code: str
    message: str


class PreflightResult(BaseModel):
    """The 200 body of ``POST /v1/actions/build/preflight``; ``ok`` says whether the build would start."""

    ok: bool
    errors: list[PreflightError] = []
    resolved: dict[str, Any] = {}


# ---------------------------------------------------------------------------
# Local files (``/v1/files``)
# ---------------------------------------------------------------------------


class LocalFileEntry(BaseModel):
    """One file or directory inside an asset's content directory."""

    #: POSIX path relative to the asset directory, e.g. ``index/hg38.1.bt2``.
    path: str
    type: Literal["file", "dir"]
    #: Bytes for a file; ``None`` for a directory.
    size: int | None = None
    #: Whether ``GET /v1/files/assets/{digest}/preview/{path}`` will show it.
    previewable: bool = False


class LocalAssetFiles(BaseModel):
    """The 200 body of ``GET /v1/files/assets/{asset_digest}``."""

    asset_digest: str
    entries: list[LocalFileEntry]
    #: True when the listing stopped at the entry cap.
    truncated: bool = False


class LocalFilePreview(BaseModel):
    """The 200 body of ``GET /v1/files/assets/{asset_digest}/preview/{path}``."""

    path: str
    size: int
    content: str
    #: True when ``content`` is only the first ``limit`` bytes of the file.
    truncated: bool
    limit: int
