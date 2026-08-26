"""The jobs wire contract: what the web UI codes against.

Everything here is a **wire** model. The internal record the manager mutates
under its lock is `Job`, which lives in `refgenie.server.jobs.manager` beside
the only code allowed to touch it; `JobManager` snapshots a `Job` into a
`JobRecord` so a reader (an SSE generator, an HTTP handler) never observes a
half-updated job.

The parameter models are thin pydantic models over the `Refgenie` facade, not
the CLI's command models: the web has different defaults (`force_large` is
always on, prompts are always refused) and a different surface (no `push_to`,
no docker flags in v1).
"""

from datetime import datetime, timezone
from enum import StrEnum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, model_validator

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
    "JobSubmission",
    "JobTarget",
    "PullJobParams",
    "describe",
    "target_of",
]


def utcnow() -> datetime:
    """Timezone-aware UTC now. One definition so every timestamp agrees."""
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
TERMINAL_STATUSES = frozenset(
    {JobStatus.SUCCEEDED, JobStatus.FAILED, JobStatus.CANCELLED}
)
#: The two from which something still might.
ACTIVE_STATUSES = frozenset({JobStatus.QUEUED, JobStatus.RUNNING})


class JobStatusFilter(StrEnum):
    """What `GET /v1/jobs?status=` accepts.

    The five real statuses, plus the two groupings a console actually asks
    for. `active` and `terminal` are not statuses and never appear on a
    record -- they exist because "what is running" and "what just finished"
    are one question each, not two and three.
    """

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
    """Base for every parameter model: unknown fields are an error.

    `extra="forbid"` is deliberate. A typo'd field name that is silently
    dropped becomes a job that quietly does the wrong thing, hours later.
    """

    model_config = ConfigDict(extra="forbid")


class PullJobParams(JobParams):
    """Parameters for a pull job. Mirrors `Refgenie.pull`'s web-relevant half."""

    server_url: str | None = None
    genome_name: str | None = None
    genome_digest: str | None = None
    asset_group_name: str
    asset_name: str | None = None
    #: Never None. A tri-state `force` means "ask the user", and there is no
    #: user to ask on a worker thread -- see `runners.pull_runner`.
    force: bool = False

    @model_validator(mode="after")
    def _genome_identified(self) -> "PullJobParams":
        if not self.genome_name and not self.genome_digest:
            raise ValueError("One of genome_name or genome_digest is required.")
        return self


class BuildParamsSpec(JobParams):
    """User-supplied build inputs, mirroring `refgenie.models.BuildParams`.

    A separate (JSON-plain) model rather than `BuildParams` itself: job
    parameters must serialize stably for `params_key` hashing and job records,
    and `BuildParams.files` holds `Path` objects. The runner converts.
    """

    assets: dict[str, str] | None = None
    params: dict[str, str | int | float | bool] | None = None
    #: Server-local paths. Acceptable only because the actions router is
    #: local-mode-only, loopback-bound, and header-guarded.
    files: dict[str, str] | None = None


class BuildJobParams(JobParams):
    """Parameters for a build job.

    `push_to` is deliberately absent: v1 has no remotes UI, and a build that
    silently registers push intent would be a surprise.
    """

    recipe_name: str
    genome_name: str
    asset_group_name: str
    asset_name: str | None = None
    recipe_version: str | None = None
    asset_description: str | None = None
    stage: bool = False
    pull_parents: bool = False
    params: BuildParamsSpec | None = None


class GenomeInitJobParams(JobParams):
    """Parameters for initializing a genome from a FASTA, then building it."""

    fasta: str
    aliases: list[str]
    description: str | None = None
    species: str | None = None
    build_fasta_asset: bool = True

    @model_validator(mode="after")
    def _at_least_one_alias(self) -> "GenomeInitJobParams":
        if not self.aliases:
            raise ValueError("At least one alias is required.")
        return self


# ---------------------------------------------------------------------------
# Wire models
# ---------------------------------------------------------------------------


#: The SSE event names. Exactly six -- this list is mirrored by
#: `JOB_EVENT_TYPES` in `frontend/src/services/contracts.ts`, and the client
#: registers one `addEventListener` per entry. An event named anything else is
#: never delivered to JS at all, so a seventh name here is a silent data loss,
#: not a compatible extension.
JobEventType = Literal["status", "log", "progress", "done", "truncated", "heartbeat"]


class JobEvent(BaseModel):
    """One entry in the manager's event log, and one SSE frame.

    One model for six frame types, so every field but `seq`/`ts`/`type` is
    optional. The router serializes with `exclude_none`, which is what keeps a
    `log` frame from carrying six null progress fields.

    `seq` is **manager-global**: there is one multiplexed stream over every
    job, so a per-job sequence could not be resumed against. Every event
    therefore carries its `job_id` -- except `truncated`, which describes the
    manager's ring rather than any one job.
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
    # These are the `JobProgress` fields, flattened. The client destructures
    # the frame and merges the remainder into the job's progress object, so
    # they must be top-level and must use these exact names.
    phase: str | None = None
    percent: float | None = None
    bytes_done: int | None = None
    bytes_total: int | None = None

    # -- log frames ------------------------------------------------------
    #: One line. NOT `message`: the client reads `lines ?? [line]` and silently
    #: drops a frame that carries neither.
    line: str | None = None
    #: "refgenie" for the package logger, "pipeline" for tee'd subprocess
    #: output. The console renders the two differently.
    source: str | None = None
    level: str | None = None

    # -- truncated frames ------------------------------------------------
    dropped: int | None = None


class JobProgress(BaseModel):
    """The most recent progress reading, denormalized onto the job record.

    A client that connects late gets the current bar position from
    `GET /v1/jobs` without replaying the whole event log.

    `phase` is a name from a shared vocabulary (see
    `frontend/src/components/jobs/phases.ts`), not free text: three of the
    phases are multi-minute silent hashing passes, and the UI shows explicit
    "this is normal and slow" copy for them. An unrecognized phase degrades to
    a humanized label rather than breaking.
    """

    #: Never null -- a progress reading with no phase has nothing to say.
    phase: str
    #: 0-100. None for every build and any pull with no Content-Length: an
    #: indeterminate bar, not a zero-length one.
    percent: float | None = None
    message: str | None = None
    bytes_done: int | None = None
    bytes_total: int | None = None


class JobTarget(BaseModel):
    """What a job acts on, structurally.

    This is the identity the UI dedupes on: a "Pull" button asks whether any
    active job already targets this asset, and a pull submitted by alias and
    one submitted by digest are the same work. Hence `genome_digest` and
    `genome_name` both, with the digest preferred.
    """

    genome_digest: str | None = None
    genome_name: str | None = None
    asset_group_name: str
    asset_name: str | None = None


class JobResult(BaseModel):
    """What a finished job produced.

    `asset_digest` and `registry_path` are what the UI renders and links to;
    the rest is there so a client can refresh exactly the rows that changed
    instead of reloading everything. Every field is nullable because
    `genome_init` produces a genome, not an asset.
    """

    asset_digest: str | None = None
    registry_path: str | None = None
    asset_group_name: str | None = None
    asset_name: str | None = None
    genome_digest: str | None = None
    staged: bool = False
    created: bool | None = None
    aliases: list[str] | None = None


class JobError(BaseModel):
    """Why a job failed, in the shared vocabulary.

    `detail` carries a traceback. That is safe here and nowhere else: local
    mode binds 127.0.0.1 and serves exactly one user, who owns the process.

    `field` names the request field at fault when there is one, so a form can
    put the message next to the input instead of in a toast.
    """

    code: ErrorCode
    message: str
    detail: str | None = None
    field: str | None = None


class JobLinks(BaseModel):
    """Where to look next. `events` is the ONE multiplexed stream, not a
    per-job one -- browsers cap HTTP/1.1 at ~6 connections per origin."""

    self: str
    events: str
    cancel: str


class JobRef(BaseModel):
    """The immediate answer to a submission. Also the 202 response body.

    Deliberately NOT the parent of `JobRecord`, and deliberately keyed
    `job_id` where the record is keyed `id`. The two are different messages:
    a ref is a receipt for a submission, a record is the state of a job. The
    client's `JobRef` and `Job` types make the same distinction, and collapsing
    them here would rename one of the two on the wire.
    """

    job_id: str
    kind: JobKind
    status: JobStatus
    created_at: datetime
    #: True when an identical job was already in flight and was returned
    #: instead of starting a second one. The UI focuses the existing card.
    duplicate: bool = False
    links: JobLinks


class JobRecord(BaseModel):
    """The full read model for one job.

    Keyed `id`, because that is what the client's job map is keyed on and what
    its "is this a whole record or a patch?" test looks for.
    """

    id: str
    kind: JobKind
    status: JobStatus
    #: Human-readable summary of the work, e.g. "pull rCRSd/fasta". The client
    #: renders a provisional one on submit and replaces it with this.
    label: str
    #: Structural identity, for the duplicate/active-job lookup.
    target: JobTarget
    created_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    #: Position in THIS job's own executor queue (0 = next up), None once
    #: running or terminal. Pulls and builds are two independent queues, so a
    #: global position would be meaningless.
    queue_position: int | None = None
    progress: JobProgress | None = None
    result: JobResult | None = None
    error: JobError | None = None
    #: Whether cancelling would actually do something. A running build is not
    #: cancellable, so the UI disables the button instead of lying.
    cancellable: bool = False
    #: How many log lines this job has produced, so the console can show a
    #: count without holding the lines.
    log_lines: int = 0
    #: True when a build succeeded without running anything, because the asset
    #: already existed. None when not applicable or not knowable.
    skipped: bool | None = None
    #: The resolved build commands, when the record carries them.
    build_commands: list[str] | None = None
    #: Not part of the client's `Job` type; kept for the API and for `curl`.
    params: dict[str, Any] = {}
    event_seq: int = 0
    log_file: str | None = None
    links: JobLinks


class JobSubmission(BaseModel):
    """`POST /v1/jobs` body: a discriminated union on `kind`.

    The actions router does not use this -- it has one typed endpoint per
    action and calls `submit_pull` / `submit_build` directly. This exists so
    the jobs API is complete and drivable by hand with curl.
    """

    kind: JobKind
    params: dict[str, Any]


#: kind -> the parameter model that validates its `params`.
PARAMS_MODEL: dict[JobKind, type[JobParams]] = {
    JobKind.PULL: PullJobParams,
    JobKind.BUILD: BuildJobParams,
    JobKind.GENOME_INIT: GenomeInitJobParams,
}


def describe(kind: JobKind, params: JobParams) -> str:
    """A one-line human label for a job.

    Mirrors the client's provisional label (`pullLabel` in
    `hooks/usePullAction.ts`) so the card does not visibly rewrite itself the
    moment the first record arrives.
    """
    if kind == JobKind.PULL:
        genome = params.genome_digest or params.genome_name or ""
        asset = f":{params.asset_name}" if params.asset_name else ""
        source = f" from {params.server_url}" if params.server_url else ""
        return f"pull {genome}/{params.asset_group_name}{asset}{source}"
    if kind == JobKind.BUILD:
        asset = f":{params.asset_name}" if params.asset_name else ""
        return f"build {params.genome_name}/{params.asset_group_name}{asset}"
    return f"initialize {params.aliases[0]}"


def target_of(kind: JobKind, params: JobParams) -> JobTarget:
    """The structural identity of the work a job does.

    `genome_init` reports the fasta asset group it ends up building, so a
    concurrent "build fasta for this genome" is recognized as the same target.
    """
    if kind == JobKind.PULL:
        return JobTarget(
            genome_digest=params.genome_digest,
            genome_name=params.genome_name,
            asset_group_name=params.asset_group_name,
            asset_name=params.asset_name,
        )
    if kind == JobKind.BUILD:
        return JobTarget(
            genome_name=params.genome_name,
            asset_group_name=params.asset_group_name,
            asset_name=params.asset_name,
        )
    return JobTarget(genome_name=params.aliases[0], asset_group_name="fasta")
