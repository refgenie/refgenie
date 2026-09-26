/**
 * The job dock: one persistent surface at the bottom of the viewport where
 * every long operation reports. The initiating click gets a one-line "queued"
 * toast and nothing else — all subsequent feedback lives here, so the user
 * never has to stay on the page that started the work.
 */

import { useEffect, useMemo, useRef, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useJobStore, selectQueued, selectRecentTerminal, selectRunning } from '../../stores/jobStore';
import { cn } from '../../utils/cn';
import { JobCard } from './JobCard';
import type { JobErrorActionId } from './errorPresentation';
import type { ConnectionState } from '../../stores/jobStore';
import type { Job } from '../../services/contracts';

const COLLAPSED_KEY = 'refgenie.job-console.collapsed';

const CONNECTION_COPY: Record<ConnectionState, { label: string; title: string }> = {
  connecting: { label: 'connecting', title: 'connecting to the live update stream' },
  open: { label: 'live', title: 'live updates' },
  polling: { label: 'polling', title: 'updating every 2s (live stream unavailable)' },
  offline: { label: 'offline', title: 'not connected' },
};

function readCollapsed(): boolean {
  try {
    return globalThis.localStorage?.getItem(COLLAPSED_KEY) === '1';
  } catch {
    return false;
  }
}

function writeCollapsed(value: boolean): void {
  try {
    globalThis.localStorage?.setItem(COLLAPSED_KEY, value ? '1' : '0');
  } catch {
    // A blocked localStorage is not a reason to break the console.
  }
}

export function JobConsole() {
  const jobs = useJobStore((state) => state.jobs);
  const connection = useJobStore((state) => state.connection);
  const clearTerminal = useJobStore((state) => state.clearTerminal);
  const navigate = useNavigate();

  const [collapsed, setCollapsed] = useState(readCollapsed);
  const [announcement, setAnnouncement] = useState('');
  const knownIds = useRef<Set<string>>(new Set());
  const knownStatuses = useRef<Map<string, string>>(new Map());

  const running = useMemo(() => selectRunning(jobs), [jobs]);
  const queued = useMemo(() => selectQueued(jobs), [jobs]);
  const recent = useMemo(() => selectRecentTerminal(jobs), [jobs]);

  // A new job auto-expands the console; it then stays expanded until the user
  // collapses it themselves.
  useEffect(() => {
    let appeared = false;
    for (const id of Object.keys(jobs)) {
      if (!knownIds.current.has(id)) {
        knownIds.current.add(id);
        appeared = true;
      }
    }
    if (appeared) {
      setCollapsed(false);
      writeCollapsed(false);
    }
  }, [jobs]);

  // Terminal transitions get a one-line spoken status. Log lines never do:
  // announcing every line would make the page unusable with a screen reader.
  useEffect(() => {
    for (const job of Object.values(jobs)) {
      const previous = knownStatuses.current.get(job.id);
      knownStatuses.current.set(job.id, job.status);
      if (previous === job.status) continue;
      if (job.status === 'succeeded') setAnnouncement(`${job.label} finished`);
      else if (job.status === 'failed') setAnnouncement(`${job.label} failed`);
      else if (job.status === 'cancelled') setAnnouncement(`${job.label} was cancelled`);
    }
  }, [jobs]);

  const visible: Job[] = [...running, ...queued, ...recent];
  const status = CONNECTION_COPY[connection];

  const handleErrorAction = (action: JobErrorActionId) => {
    // The console is not the place these are fixed; it hands off to the screen
    // that owns the form, rather than growing a second copy of it.
    if (action === 'fix_inputs' || action === 'build_input') navigate('/build');
    else if (action === 'init_genome' || action === 'manage_recipes' || action === 'pick_server')
      navigate('/manage');
    else navigate('/remote');
  };

  if (visible.length === 0 && connection === 'open') return null;

  const toggle = () => {
    setCollapsed((current) => {
      writeCollapsed(!current);
      return !current;
    });
  };

  return (
    <aside
      className={cn('rg-job-console', collapsed && 'rg-job-console--collapsed')}
      aria-label="Job console"
    >
      <header className="rg-job-console__header">
        <button
          type="button"
          className="rg-job-console__toggle"
          onClick={toggle}
          aria-expanded={!collapsed}
        >
          <span aria-hidden="true">{collapsed ? '▲' : '▼'}</span> Jobs
        </button>

        <span className="rg-job-console__badge" aria-live="polite">
          {running.length} running, {queued.length} queued
        </span>

        <span
          className={cn('rg-job-console__status', `rg-job-console__status--${connection}`)}
          title={status.title}
        >
          <span className="rg-job-console__dot" aria-hidden="true" />
          {status.label}
        </span>

        <span className="sr-only" aria-live="polite">
          {announcement}
        </span>

        {recent.length > 0 && (
          <button type="button" className="rg-btn rg-btn--sm rg-btn--bare" onClick={clearTerminal}>
            Clear finished
          </button>
        )}
      </header>

      {!collapsed && (
        <div className="rg-job-console__body">
          {visible.length === 0 ? (
            <p className="rg-muted text-sm">No jobs running.</p>
          ) : (
            visible.map((job) => (
              <JobCard key={job.id} job={job} onErrorAction={handleErrorAction} />
            ))
          )}
        </div>
      )}
    </aside>
  );
}
