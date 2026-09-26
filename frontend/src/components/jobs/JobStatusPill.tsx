import { cn } from '../../utils/cn';
import type { JobStatus } from '../../services/contracts';

const LABELS: Record<JobStatus, string> = {
  queued: 'Queued',
  running: 'Running',
  succeeded: 'Succeeded',
  failed: 'Failed',
  cancelled: 'Cancelled',
};

export function JobStatusPill({ status }: { status: JobStatus }) {
  return (
    <span className={cn('rg-status-pill', `rg-status-pill--${status}`)}>{LABELS[status]}</span>
  );
}
