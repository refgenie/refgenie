/**
 * Everything about one job: the full log (from `/v1/jobs/{id}/log`, not the
 * store's 500-line live ring), the error block, the resolved build commands
 * when the record carries them, and a link to the asset the job produced.
 */

import { useEffect, useRef } from 'react';
import { Link } from 'react-router-dom';
import { BaseModal } from '../common/BaseModal';
import { DescriptionList } from '../common/DescriptionList';
import { JobErrorBlock } from './JobErrorBlock';
import { JobStatusPill } from './JobStatusPill';
import { useJobLog } from '../../hooks/queries/useJobs';
import { formatDuration, formatTimestamp } from '../../utils/time';
import type { Job } from '../../services/contracts';

export interface JobDetailModalProps {
  job: Job | null;
  isOpen: boolean;
  onClose: () => void;
  /** Scroll the log to the end on open — how a `build_failed` card arrives. */
  scrollToEnd?: boolean;
}

export function JobDetailModal({ job, isOpen, onClose, scrollToEnd }: JobDetailModalProps) {
  const log = useJobLog(job?.id, 0, { enabled: isOpen && !!job });
  const logRef = useRef<HTMLPreElement>(null);

  const lineCount = log.data?.lines.length ?? 0;
  useEffect(() => {
    if (!isOpen || !scrollToEnd) return;
    const node = logRef.current;
    if (node) node.scrollTop = node.scrollHeight;
  }, [isOpen, scrollToEnd, lineCount]);

  if (!job) return null;

  return (
    <BaseModal isOpen={isOpen} onClose={onClose} title={job.label} size="lg">
      <div className="flex flex-col gap-4">
        <DescriptionList
          items={[
            { term: 'Status', value: <JobStatusPill status={job.status} /> },
            { term: 'Kind', value: job.kind },
            { term: 'Job id', value: <code className="rg-code rg-code--inline">{job.id}</code> },
            { term: 'Created', value: formatTimestamp(job.created_at) },
            { term: 'Started', value: formatTimestamp(job.started_at) },
            { term: 'Finished', value: formatTimestamp(job.finished_at) },
            { term: 'Duration', value: formatDuration(job.started_at, job.finished_at) },
            ...(job.result
              ? [
                  {
                    term: 'Result',
                    value: (
                      <Link className="rg-link" to={`/assets/${job.result.asset_digest}`}>
                        {job.result.registry_path}
                      </Link>
                    ),
                  },
                ]
              : []),
          ]}
        />

        <JobErrorBlock job={job} />

        {job.build_commands && job.build_commands.length > 0 && (
          <section>
            <h3 className="text-sm font-semibold mb-2">Build commands</h3>
            <pre className="rg-code rg-job-detail__log">{job.build_commands.join('\n')}</pre>
          </section>
        )}

        <section>
          <h3 className="text-sm font-semibold mb-2">Log</h3>
          {log.isPending ? (
            <p className="rg-muted text-sm">Loading log…</p>
          ) : log.error ? (
            <p className="rg-muted text-sm">
              No log is available for this job. Jobs live in the server process, so a
              restart takes their logs with them.
            </p>
          ) : lineCount === 0 ? (
            <p className="rg-muted text-sm">This job produced no log output.</p>
          ) : (
            <pre className="rg-code rg-job-detail__log" ref={logRef}>
              {log.data?.truncated ? '… earlier lines dropped …\n' : ''}
              {log.data?.lines.join('\n')}
            </pre>
          )}
        </section>

        <BaseModal.Footer>
          <button type="button" className="rg-btn" onClick={onClose}>
            Close
          </button>
        </BaseModal.Footer>
      </div>
    </BaseModal>
  );
}
