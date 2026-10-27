"""In-process background jobs for the local refgenie web UI.

A pull takes minutes and a build takes hours; neither can happen inside a
request handler. This package runs them on dedicated worker threads, records
status, progress, log lines, result and error per job, and streams that state
to the browser over one multiplexed SSE endpoint (with a JSON polling
fallback).

The names re-exported here are the public surface. The sibling actions router
imports from this package rather than reaching into its submodules::

    from refgenie.server.jobs import (
        JobManager, JobRef, PullJobParams, get_job_manager,
    )

    @router.post("/pull", status_code=202, response_model=JobRef)
    def pull(req: PullRequest, jobs: JobManager = Depends(get_job_manager)) -> JobRef:
        return jobs.submit_pull(PullJobParams(genome_name=req.genome, ...))

`submit_pull` / `submit_build` / `submit_genome_init` never block and never
raise for a job-level failure -- a failed pull is a job record with an
`error.code`, not an HTTP 500.
"""

from refgenie.server.jobs.manager import CancelOutcome, JobContext, JobManager
from refgenie.server.jobs.schemas import (
    BuildJobParams,
    GenomeInitJobParams,
    JobError,
    JobEvent,
    JobKind,
    JobRecord,
    JobRef,
    JobStatus,
    PullJobParams,
)
from refgenie.server.jobs.router import get_job_manager
from refgenie.server.jobs.router import router as jobs_router

__all__ = [
    "BuildJobParams",
    "CancelOutcome",
    "GenomeInitJobParams",
    "JobContext",
    "JobError",
    "JobEvent",
    "JobKind",
    "JobManager",
    "JobRecord",
    "JobRef",
    "JobStatus",
    "PullJobParams",
    "get_job_manager",
    "jobs_router",
]
