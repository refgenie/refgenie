"""Web request models for the local actions API, and its response envelopes.

These are **web** models, not the CLI's command models: the CLI models key
their JSON schemas by CLI aliases (one-letter flags included) and mostly lack
``populate_by_name``, so a field-name JSON body fails validation outright. The
web contract is snake_case field names, no aliases, no ``validation_alias``.

``extra="forbid"`` on every request model is deliberate: it turns a typo'd or
stale field name into a 422 instead of a silent no-op, which is exactly the
class of rot that hid the broken dash pages.

Two deliberate design points:

* ``AliasSetRequest.genome_digest`` is **required**. The facade's digest-less
  variant performs synchronous network lookups against every subscription; a
  request/response endpoint must not do that. The UI always knows the digest.
* ``BuildParamsRequest.files`` and ``GenomeInitRequest.fasta`` accept
  server-local paths. That is a real capability the build path needs, and it
  is acceptable only because local mode binds 127.0.0.1, requires the action
  header, and enforces the host guard. Never expose this router in server
  mode.

This module imports only pydantic -- no manager imports -- so it stays cheap
and testable.
"""

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

__all__ = [
    "ActionRequest",
    "ActionResult",
    "AliasSetRequest",
    "BuildParamsRequest",
    "BuildRequest",
    "GenomeInitRequest",
    "PreflightError",
    "PreflightResult",
    "PullRequest",
    "SetDefaultAssetRequest",
    "SubscribeRequest",
    "UnsubscribeRequest",
]


class ActionRequest(BaseModel):
    """Base for every request model: unknown fields are an error."""

    model_config = ConfigDict(extra="forbid")


class PullRequest(ActionRequest):
    """Body of ``POST /v1/actions/pull``."""

    asset_group: str
    genome: str | None = None  # alias name
    genome_digest: str | None = None
    asset: str | None = None  # specific asset name within the group
    server_url: str | None = None  # restrict the pull to one server
    #: Never None: a tri-state force means "ask the user", and the web layer
    #: never prompts.
    force: bool = False

    @model_validator(mode="after")
    def _exactly_one_genome_ref(self) -> "PullRequest":
        if bool(self.genome) == bool(self.genome_digest):
            raise ValueError("Exactly one of 'genome' or 'genome_digest' must be set.")
        return self


class BuildParamsRequest(ActionRequest):
    """User-supplied build inputs (see the module docstring about ``files``)."""

    assets: dict[str, str] | None = None
    params: dict[str, str | int | float | bool] | None = None
    files: dict[str, str] | None = None


class BuildRequest(ActionRequest):
    """Body of ``POST /v1/actions/build`` and ``/build/preflight``."""

    recipe: str
    genome: str
    asset_group: str
    asset: str | None = None
    recipe_version: str | None = None
    description: str | None = None
    stage: bool = False
    pull_parents: bool = False
    params: BuildParamsRequest | None = None


class AliasSetRequest(ActionRequest):
    """Body of ``POST /v1/actions/aliases``."""

    alias: str
    genome_digest: str


class SubscribeRequest(ActionRequest):
    """Body of ``POST /v1/actions/subscriptions``."""

    server_urls: list[str] = Field(min_length=1)
    reset: bool = False


class UnsubscribeRequest(ActionRequest):
    """Body of ``DELETE /v1/actions/subscriptions``."""

    server_urls: list[str] = Field(min_length=1)


class SetDefaultAssetRequest(ActionRequest):
    """Body of ``POST /v1/actions/assets/default``."""

    genome_digest: str
    asset_group: str
    asset: str


class GenomeInitRequest(ActionRequest):
    """Body of ``POST /v1/actions/genomes``."""

    fasta: str  # server-local path or URL
    aliases: list[str] = Field(min_length=1)
    description: str | None = None
    species: str | None = None
    build_fasta_asset: bool = True  # ingest, then build the fasta asset


# ---------------------------------------------------------------------------
# Response envelopes
# ---------------------------------------------------------------------------
#
# There is no JobAccepted envelope: a 202 returns the jobs package's `JobRef`
# directly ({job_id, kind, status, created_at, duplicate, links}).
# `duplicate: true` means an identical in-flight job was coalesced, not an
# error. Non-2xx bodies are `refgenie.server.errors.ErrorResponse`.


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
    """The 200 body of ``POST /v1/actions/build/preflight``.

    Always 200: a preflight that *found problems* has succeeded at its job.
    ``ok`` reports whether the build would start.
    """

    ok: bool
    errors: list[PreflightError] = []
    resolved: dict[str, Any] = {}
