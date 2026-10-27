/**
 * The failure surface of a job card, driven entirely by `error.code`.
 *
 * The offered action is an id, not a callback baked into the table: the same
 * `pull_force` presentation means "re-submit with force" on a pull card and
 * "not available" on the history page, and the owner decides.
 */

import { useJobStore } from '../../stores/jobStore';
import { CopyButton } from '../common/CopyButton';
import { cn } from '../../utils/cn';
import { presentJobError } from './errorPresentation';
import type { JobErrorActionId } from './errorPresentation';
import type { Job } from '../../services/contracts';

export interface JobErrorBlockProps {
  job: Job;
  /** Return false (or omit the prop) to hide the offered action. */
  onAction?: (action: JobErrorActionId, job: Job) => void;
  availableActions?: readonly JobErrorActionId[];
}

/** The blob a user pastes into a GitHub issue. */
function buildReport(job: Job, logLines: string[]): string {
  return [
    `job: ${job.id}`,
    `kind: ${job.kind}`,
    `status: ${job.status}`,
    `target: ${JSON.stringify(job.target)}`,
    `code: ${job.error?.code ?? 'none'}`,
    `message: ${job.error?.message ?? 'none'}`,
    `detail: ${job.error?.detail ?? 'none'}`,
    '',
    'last log lines:',
    ...logLines,
  ].join('\n');
}

export function JobErrorBlock({ job, onAction, availableActions }: JobErrorBlockProps) {
  const logs = useJobStore((state) => state.logs[job.id]);
  if (!job.error) return null;

  const presentation = presentJobError(job.error);
  const tail = (logs ?? []).slice(-50).map((line) => line.text);
  const action = presentation.action;
  const actionAllowed =
    !!action && !!onAction && (!availableActions || availableActions.includes(action.id));

  return (
    <div
      className={cn('rg-job-error', `rg-job-error--${presentation.tone}`)}
      role={presentation.tone === 'error' ? 'alert' : undefined}
    >
      <p className="rg-job-error__headline">{presentation.headline}</p>

      {presentation.note && <p className="rg-job-error__note">{presentation.note}</p>}

      {presentation.preformatted ? (
        job.error.detail ? (
          <pre className="rg-code rg-job-error__log">{job.error.detail}</pre>
        ) : (
          <p className="rg-muted">{job.error.message}</p>
        )
      ) : (
        <p className="rg-job-error__message">{job.error.message}</p>
      )}

      {presentation.disclosure && job.error.detail && (
        <details className="rg-job-error__details">
          <summary>Details</summary>
          <pre className="rg-code rg-job-error__log">{job.error.detail}</pre>
        </details>
      )}

      <div className="rg-job-error__actions">
        {actionAllowed && action && (
          <button
            type="button"
            className="rg-btn rg-btn--sm"
            onClick={() => onAction?.(action.id, job)}
          >
            {action.label}
          </button>
        )}
        <CopyButton value={buildReport(job, tail)} label="Copy details" />
      </div>
    </div>
  );
}
