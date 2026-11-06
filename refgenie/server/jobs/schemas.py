"""The jobs wire models: parameters, records and SSE events."""

from datetime import datetime, timezone
from enum import StrEnum
from typing import Any, Literal

from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, model_validator

from refgenie.models import BuildParams
from refgenie.server.errors import ErrorCode

__all__ = [
    "BuildJobParams",
    "BuildParamsSpec",
    "GenomeInitJobParams",
    "JobError",
    "JobEvent",
    "JobEventType",
    "JobKind",
    "JobLinks",
    "JobProgress",
    "JobRecord",
    "JobRef",
    "JobResult",
    "JobStatus",
    "JobStatusFilter",
    "JobTarget",
    "PullJobParams",
    "describe",
    "target_of",
]


def utcnow() -> datetime:
    """Timezone-aware UTC now."""
    return datetime.now(timezone.utc)


class JobKind(StrEnum):
    """What a job does. One runner per kind."""

    PULL = "pull"
    BUILD = "build"
    GENOME_INIT = "genome_init"


class JobStatus(StrEnum):
    """The job state machine. Only the manager writes these."""

    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"


#: The three states from which nothing further happens.
TERMINAL_STATUSES = frozenset({JobStatus.SUCCEEDED, JobStatus.FAILED, JobStatus.CANCELLED})
#: The two from which something still might.
ACTIVE_STATUSES = frozenset({JobStatus.QUEUED, JobStatus.RUNNING})


class JobStatusFilter(StrEnum):
    """What `GET /v1/jobs?status=` accepts: a status, `active` or `terminal`."""

    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"
    ACTIVE = "active"
    TERMINAL = "terminal"

    def statuses(self) -> frozenset[JobStatus]:
        """The concrete statuses this filter selects."""
        if self is JobStatusFilter.ACTIVE:
            return ACTIVE_STATUSES
        if self is JobStatusFilter.TERMINAL:
            return TERMINAL_STATUSES
        return frozenset({JobStatus(self.value)})


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------


class JobParams(BaseModel):
    """Base for every parameter model: unknown fields are an error."""

    model_config = ConfigDict(extra="forbid")


class PullJobParams(JobParams):
    """Body of `POST /v1/actions/pull`, and a pull job's parameters."""

    asset_group: str
    #: An alias; set exactly one of `genome` or `genome_digest`.
    genome: str | None = None
    genome_digest: str | None = None
    asset: str | None = None
    #: Restrict the pull to one server.
    server_url: str | None = None
    #: Never None: a worker thread has no user to ask.
    force: bool = False

    @model_validator(mode="after")
    def _exactly_one_genome_ref(self) -> "PullJobParams":
        if bool(self.genome) == bool(self.genome_digest):
            raise ValueError("Exactly one of 'genome' or 'genome_digest' must be set.")
        return self


class BuildParamsSpec(JobParams):
    """User-supplied build inputs, JSON-plain."""

    assets: dict[str, str] | None = None
    params: dict[str, str | int | float | bool] | None = None
    #: Server-local paths; local mode only.
    files: dict[str, str] | None = None

    def to_build_params(self) -> BuildParams:
        """The library-side `BuildParams`."""
        files = {name: Path(path) for name, path in self.files.items()} if self.files else None
        return BuildParams(assets=self.assets, params=self.params, files=files)


class BuildJobParams(JobParams):
    """Body of `POST /v1/actions/build` and `/build/preflight`, and a build job's parameters."""

    recipe: str
    #: An alias; the build folder is named after it.
    genome: str
    asset_group: str
    asset: str | None = None
    recipe_version: str | None = None
    description: str | None = None
    stage: bool = False
    pull_parents: bool = False
    params: BuildParamsSpec | None = None


class GenomeInitJobParams(JobParams):
    """Body of `POST /v1/actions/genomes`: load a FASTA, then build its fasta asset."""

    #: A server-local path or URL.
    fasta: str
    aliases: list[str] = Field(min_length=1)
    description: str | None = None
    species: str | None = None
    build_fasta_asset: bool = True


# ---------------------------------------------------------------------------
# Wire models
# ---------------------------------------------------------------------------


#: The SSE event names, mirrored by `JOB_EVENT_TYPES` in
#: `frontend/src/services/contracts.ts`. A new name is not a compatible change.
JobEventType = Literal["status", "log", "progress", "done", "truncated", "heartbeat"]


class JobEvent(BaseModel):
    """One entry in the manager's event log, and one SSE frame.

    `seq` is manager-global. Every event but `truncated` carries a `job_id`.
    """

    seq: int
    ts: datetime
    type: JobEventType
    job_id: str | None = None

    # -- status frames ---------------------------------------------------
    status: JobStatus | None = None
    queue_position: int | None = None
    message: str | None = None

    # -- progress frames -------------------------------------------------
    # The `JobProgress` fields, flattened; the client merges them by name.
    phase: str | None = None
    percent: float | None = None
    bytes_done: int | None = None
    bytes_total: int | None = None

    # -- log frames ------------------------------------------------------
    #: One log line (the client reads `line`, not `message`).
    line: str | None = None
    #: "refgenie" for the package logger, "pipeline" for subprocess output.
    source: str | None = None
    level: str | None = None

    # -- truncated frames ------------------------------------------------
    dropped: int | None = None


class JobProgress(BaseModel):
    """The most recent progress reading, kept on the job record.

    `phase` is a name from `frontend/src/components/jobs/phases.ts`.
    """

    phase: str
    #: 0-100, or None when the total is unknown.
    percent: float | None = None
    message: str | None = None
    bytes_done: int | None = None
    bytes_total: int | None = None


class JobTarget(BaseModel):
    """What a job acts on; the UI matches active jobs to assets by this."""

    genome_digest: str | None = None
    genome_name: str | None = None
    asset_group_name: str
    asset_name: str | None = None


class JobResult(BaseModel):
    """What a finished job produced; asset fields are empty for `genome_init`."""

    asset_digest: str | None = None
    registry_path: str | None = None
    asset_group_name: str | None = None
    asset_name: str | None = None
    genome_digest: str | None = None
    staged: bool = False
    created: bool | None = None
    aliases: list[str] | None = None


class JobError(BaseModel):
    """Why a job failed.

    `detail` is a traceback (local mode only). `field` names the request field
    at fault, when there is one.
    """

    code: ErrorCode
    message: str
    detail: str | None = None
    field: str | None = None


class JobLinks(BaseModel):
    """Links for a job. `events` is the shared multiplexed stream."""

    self: str
    events: str
    cancel: str


class JobRef(BaseModel):
    """The receipt for a submission, keyed `job_id` (a `JobRecord` is keyed `id`)."""

    job_id: str
    kind: JobKind
    status: JobStatus
    created_at: datetime
    #: True when an identical active job was returned instead of a new one.
    duplicate: bool = False
    links: JobLinks


class JobRecord(BaseModel):
    """The full read model for one job."""

    id: str
    kind: JobKind
    status: JobStatus
    #: Human-readable summary of the work, e.g. "pull rCRSd/fasta".
    label: str
    target: JobTarget
    created_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    #: Position in this job's executor queue (0 = next up); None once running.
    queue_position: int | None = None
    progress: JobProgress | None = None
    result: JobResult | None = None
    error: JobError | None = None
    #: Whether cancelling would do anything (false for a running build).
    cancellable: bool = False
    #: How many log lines this job has produced.
    log_lines: int = 0
    #: True when a build succeeded without running because the asset existed.
    skipped: bool | None = None
    params: dict[str, Any] = {}
    event_seq: int = 0
    log_file: str | None = None
    links: JobLinks


def describe(kind: JobKind, params: JobParams) -> str:
    """A one-line human label for a job, matching the client's provisional label."""
    if kind == JobKind.PULL:
        genome = params.genome_digest or params.genome
        asset = f":{params.asset}" if params.asset else ""
        source = f" from {params.server_url}" if params.server_url else ""
        return f"pull {genome}/{params.asset_group}{asset}{source}"
    if kind == JobKind.BUILD:
        asset = f":{params.asset}" if params.asset else ""
        return f"build {params.genome}/{params.asset_group}{asset}"
    return f"initialize {params.aliases[0]}"


def target_of(kind: JobKind, params: JobParams) -> JobTarget:
    """The structural identity of the work a job does (`genome_init` targets fasta)."""
    if kind == JobKind.PULL:
        return JobTarget(
            genome_digest=params.genome_digest,
            genome_name=params.genome,
            asset_group_name=params.asset_group,
            asset_name=params.asset,
        )
    if kind == JobKind.BUILD:
        return JobTarget(
            genome_name=params.genome,
            asset_group_name=params.asset_group,
            asset_name=params.asset,
        )
    return JobTarget(genome_name=params.aliases[0], asset_group_name="fasta")
