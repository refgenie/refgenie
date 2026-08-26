/**
 * One job.
 *
 * A queued card is deliberately NOT a running card with an empty bar. Builds
 * run one at a time and pulls two at a time, so waiting is the normal state,
 * and a card that showed a stalled progress bar would read as a hang. Queued
 * cards show their position and offer cancel (cancelling something that has
 * not started is always safe); running cards show progress and the log tail.
 */

import { useEffect, useState } from 'react';
import { useCancelJob } from '../../hooks/queries/useJobs';
import { useCapability } from '../../hooks/useCapability';
import { useJobStore } from '../../stores/jobStore';
import { cn } from '../../utils/cn';
import { formatDuration } from '../../utils/time';
import { isTerminal } from '../../services/contracts';
import { JobProgress } from './JobProgress';
import { JobLogTail } from './JobLogTail';
import { JobErrorBlock } from './JobErrorBlock';
import { JobStatusPill } from './JobStatusPill';
import { JobDetailModal } from './JobDetailModal';
import type { JobErrorActionId } from './errorPresentation';
import type { Job } from '../../services/contracts';

const KIND_ICON: Record<Job['kind'], string> = {
  pull: '↓',
  build: '⚙',
  genome_init: '+',
};

const KIND_LABEL: Record<Job['kind'], string> = {
  pull: 'Pull',
  build: 'Build',
  genome_init: 'Initialize genome',
};

/** Builds have no percentage at all, so their log tail carries the whole story. */
const LOG_LINES: Record<Job['kind'], number> = { pull: 12, build: 20, genome_init: 12 };

export interface JobCardProps {
  job: Job;
  onErrorAction?: (action: JobErrorActionId, job: Job) => void;
}

function queuedLabel(position: number | null): string {
  if (position === null || position === undefined) return 'Queued';
  if (position <= 0) return 'Queued — next up';
  return `Queued — ${position} ahead`;
}

export function JobCard({ job, onErrorAction }: JobCardProps) {
  const canCancel = useCapability('jobs_cancel');
  const cancel = useCancelJob();
  const focusedJobId = useJobStore((state) => state.focusedJobId);
  const [detailOpen, setDetailOpen] = useState(false);
  const [now, setNow] = useState(() => Date.now());

  const running = job.status === 'running';
  useEffect(() => {
    if (!running) return;
    const timer = setInterval(() => setNow(Date.now()), 1000);
    return () => clearInterval(timer);
  }, [running]);

  const queued = job.status === 'queued';
  const terminal = isTerminal(job.status);
  // A build whose asset already existed succeeds without doing any work. Saying
  // "succeeded" alone would imply something was rebuilt.
  const skipped = job.status === 'succeeded' && job.skipped === true;

  const handleErrorAction = (action: JobErrorActionId, target: Job) => {
    if (action === 'view_log') {
      setDetailOpen(true);
      return;
    }
    onErrorAction?.(action, target);
  };

  return (
    <article
      className={cn(
        'rg-job-card',
        queued && 'rg-job-card--queued',
        job.status === 'failed' && 'rg-job-card--failed',
        job.id === focusedJobId && 'rg-job-card--focused',
      )}
    >
      <header className="rg-job-card__header">
        <span className="rg-job-card__icon" aria-hidden="true">
          {KIND_ICON[job.kind]}
        </span>
        <span className="rg-job-card__label">
          <span className="sr-only">{KIND_LABEL[job.kind]}: </span>
          {job.label}
        </span>
        <JobStatusPill status={job.status} />
        {(running || terminal) && (
          <span className="rg-job-card__elapsed rg-muted">
            {formatDuration(job.started_at, job.finished_at, now)}
          </span>
        )}
      </header>

      {queued ? (
        <p className="rg-job-card__queued">{queuedLabel(job.queue_position)}</p>
      ) : (
        <>
          {!terminal && <JobProgress kind={job.kind} progress={job.progress} />}
          {!terminal && <JobLogTail jobId={job.id} status={job.status} lines={LOG_LINES[job.kind]} />}
        </>
      )}

      {skipped && (
        <p className="rg-job-card__note">Already built — nothing to do.</p>
      )}

      {job.error && (
        <JobErrorBlock
          job={job}
          onAction={handleErrorAction}
          availableActions={
            onErrorAction
              ? undefined
              : (['view_log'] as const)
          }
        />
      )}

      <footer className="rg-job-card__actions">
        {canCancel && !terminal && (
          <button
            type="button"
            className="rg-btn rg-btn--sm"
            // A running build is not cancellable: pypiper owns the subprocess
            // tree. The record says so, and the button disables rather than lies.
            disabled={cancel.isPending || (running && job.cancellable === false)}
            title={
              running && job.cancellable === false
                ? 'A running build cannot be cancelled'
                : undefined
            }
            onClick={() => cancel.mutate(job.id)}
          >
            Cancel
          </button>
        )}
        <button
          type="button"
          className="rg-btn rg-btn--sm rg-btn--bare"
          onClick={() => setDetailOpen(true)}
        >
          Details
        </button>
      </footer>

      <JobDetailModal
        job={job}
        isOpen={detailOpen}
        onClose={() => setDetailOpen(false)}
        scrollToEnd={job.error?.code === 'build_failed'}
      />
    </article>
  );
}
