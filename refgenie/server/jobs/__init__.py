"""In-process background jobs (pull, build, genome init) for the local web UI.

Import from this package, not its submodules. Submitting never blocks; a failed
job is a record with an `error.code`, not an HTTP error.
"""

from refgenie.server.jobs.manager import CancelOutcome, JobContext, JobManager
from refgenie.server.jobs.router import get_job_manager
from refgenie.server.jobs.router import router as jobs_router
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
