/**
 * Job history.
 *
 * Jobs live in the server process and are intentionally not persisted, so an
 * empty history after a restart is expected rather than a bug. The banner says
 * that outright instead of leaving the user to wonder.
 */

import { useState } from 'react';
import { Link } from 'react-router-dom';
import { useCapability } from '../hooks/useCapability';
import { useJobsList } from '../hooks/queries/useJobs';
import { DataTable } from '../components/common/DataTable';
import { Pagination } from '../components/common/Pagination';
import { EmptyState } from '../components/common/states';
import { JobStatusPill } from '../components/jobs/JobStatusPill';
import { JobDetailModal } from '../components/jobs/JobDetailModal';
import { formatDuration, formatTimestamp } from '../utils/time';
import { NotAvailablePage } from './NotAvailablePage';
import type { Column } from '../components/common/DataTable';
import type { Job, JobKind } from '../services/contracts';
import type { JobStatusFilter } from '../services/jobs';

const PAGE_SIZE = 25;

const STATUS_FILTERS: Array<{ value: '' | JobStatusFilter; label: string }> = [
  { value: '', label: 'All statuses' },
  { value: 'active', label: 'Active' },
  { value: 'terminal', label: 'Finished' },
  { value: 'queued', label: 'Queued' },
  { value: 'running', label: 'Running' },
  { value: 'succeeded', label: 'Succeeded' },
  { value: 'failed', label: 'Failed' },
  { value: 'cancelled', label: 'Cancelled' },
];

const KIND_FILTERS: Array<{ value: '' | JobKind; label: string }> = [
  { value: '', label: 'All kinds' },
  { value: 'pull', label: 'Pull' },
  { value: 'build', label: 'Build' },
  { value: 'genome_init', label: 'Genome init' },
];

export function JobsPage() {
  const enabled = useCapability('jobs');
  const [status, setStatus] = useState<'' | JobStatusFilter>('');
  const [kind, setKind] = useState<'' | JobKind>('');
  const [offset, setOffset] = useState(0);
  const [selected, setSelected] = useState<Job | null>(null);

  const jobs = useJobsList(
    {
      status: status || undefined,
      kind: kind || undefined,
      offset,
      limit: PAGE_SIZE,
    },
    { enabled },
  );

  if (!enabled) {
    return (
      <NotAvailablePage
        title="Jobs"
        reason="This instance does not run jobs (capability `jobs` is off)."
      />
    );
  }

  const columns: Array<Column<Job>> = [
    { key: 'status', header: 'Status', render: (job) => <JobStatusPill status={job.status} /> },
    { key: 'kind', header: 'Kind', render: (job) => job.kind },
    {
      key: 'label',
      header: 'Job',
      render: (job) => (
        <button
          type="button"
          className="rg-btn rg-btn--bare"
          onClick={() => setSelected(job)}
        >
          {job.label}
        </button>
      ),
    },
    { key: 'started', header: 'Started', render: (job) => formatTimestamp(job.started_at) },
    {
      key: 'duration',
      header: 'Duration',
      align: 'right',
      render: (job) => formatDuration(job.started_at, job.finished_at),
    },
    {
      key: 'result',
      header: 'Result',
      render: (job) =>
        job.result ? (
          <Link className="rg-link" to={`/assets/${job.result.asset_digest}`}>
            {job.result.registry_path}
          </Link>
        ) : (
          <span className="rg-muted">—</span>
        ),
    },
  ];

  return (
    <div className="flex flex-col gap-6">
      <h1 className="text-3xl font-bold">Jobs</h1>

      <div className="rg-banner rg-banner--info" role="status">
        <span>
          Jobs live in the refgenie server process and are not saved to disk. Restarting the
          server clears this history — an empty list after a restart is expected.
        </span>
      </div>

      <div className="flex gap-3 flex-wrap items-end">
        <div className="rg-field">
          <label className="rg-field__label" htmlFor="jobs-status">
            Status
          </label>
          <select
            id="jobs-status"
            className="rg-field__input"
            value={status}
            onChange={(event) => {
              setStatus(event.target.value as '' | JobStatusFilter);
              setOffset(0);
            }}
          >
            {STATUS_FILTERS.map((option) => (
              <option key={option.value} value={option.value}>
                {option.label}
              </option>
            ))}
          </select>
        </div>

        <div className="rg-field">
          <label className="rg-field__label" htmlFor="jobs-kind">
            Kind
          </label>
          <select
            id="jobs-kind"
            className="rg-field__input"
            value={kind}
            onChange={(event) => {
              setKind(event.target.value as '' | JobKind);
              setOffset(0);
            }}
          >
            {KIND_FILTERS.map((option) => (
              <option key={option.value} value={option.value}>
                {option.label}
              </option>
            ))}
          </select>
        </div>
      </div>

      <DataTable
        caption="Jobs"
        columns={columns}
        rows={jobs.data?.items}
        rowKey={(job) => job.id}
        loading={jobs.isPending}
        error={jobs.error}
        onRetry={() => jobs.refetch()}
        empty={<EmptyState message="No jobs have run since this server started." />}
      />

      <Pagination pagination={jobs.data?.pagination} onOffsetChange={setOffset} />

      <JobDetailModal
        job={selected}
        isOpen={selected !== null}
        onClose={() => setSelected(null)}
        scrollToEnd={selected?.error?.code === 'build_failed'}
      />
    </div>
  );
}
