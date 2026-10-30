/**
 * The public page's dock for pulls it handed to a connected local refgenie.
 *
 * It shares the look of the local dash's `JobConsole` (same dock classes, the
 * same status pill and progress bar) but not its card: `JobCard` reads its log
 * tail from `jobStore` and offers cancel and a log modal against THIS page's
 * backend, none of which apply to a job running on another machine. This card
 * shows progress, the outcome, and, once done, a link to the asset under
 * `/local`.
 *
 * Mounted once by AppLayout (never on a local dash). It also owns the poller,
 * so polling exists exactly as long as the dock does.
 */

import { useEffect, useRef, useState } from 'react';
import { Link } from 'react-router-dom';
import { ExternalLink } from '../common/ExternalLink';
import { Icon } from '../common/Icon';
import { DEFAULT_SERVER_URL, localServersUrl } from '../servers/links';
import { JobProgress } from '../jobs/JobProgress';
import { JobStatusPill } from '../jobs/JobStatusPill';
import { cn } from '../../utils/cn';
import { useBridge } from '../../hooks/useBridge';
import { useBridgeJobEntries, useBridgeJobPolling } from '../../hooks/useBridgeJobPolling';
import { isActiveEntry, useBridgeJobStore } from '../../stores/bridgeJobStore';
import type { BridgeJobEntry } from '../../stores/bridgeJobStore';

/** Failures fixed by subscribing, which only the local dash itself may do. */
const SUBSCRIPTION_ERRORS = new Set(['no_subscriptions', 'server_not_subscribed']);

function BridgeJobCard({ entry, localBaseUrl }: { entry: BridgeJobEntry; localBaseUrl: string | null }) {
  const dismiss = useBridgeJobStore((state) => state.dismiss);
  const job = entry.job;
  const status = job?.status ?? 'queued';
  const active = isActiveEntry(entry);

  return (
    <article
      className={cn(
        'rg-job-card',
        status === 'queued' && !entry.lost && 'rg-job-card--queued',
        status === 'failed' && 'rg-job-card--failed',
      )}
    >
      <header className="rg-job-card__header">
        <Icon name="download" className="rg-icon--lg" />
        <span className="rg-job-card__label">
          <span className="sr-only">Pull: </span>
          {entry.label}
        </span>
        {!entry.lost && <JobStatusPill status={status} />}
      </header>

      {entry.lost ? (
        <p className="rg-job-card__note">
          Your local refgenie no longer has a record of this pull, probably because it
          restarted. Check its genomes to see whether the asset arrived.
        </p>
      ) : status === 'queued' ? (
        <p className="rg-job-card__queued">Queued on your local refgenie</p>
      ) : status === 'running' ? (
        <JobProgress kind="pull" progress={job?.progress ?? null} />
      ) : status === 'succeeded' ? (
        <p className="rg-job-card__note">
          Now on your local refgenie.{' '}
          {job?.result && (
            <Link className="rg-link" to={`/local/assets/${job.result.asset_digest}`}>
              View {job.result.registry_path}
            </Link>
          )}
        </p>
      ) : status === 'failed' ? (
        <p className="rg-job-card__note">
          {job?.error?.message ?? 'The pull failed.'}
          {localBaseUrl && job?.error?.code && SUBSCRIPTION_ERRORS.has(job.error.code) && (
            <>
              {' '}
              <ExternalLink
                href={localServersUrl(localBaseUrl, entry.serverUrl ?? DEFAULT_SERVER_URL)}
              >
                Subscribe on your local refgenie
              </ExternalLink>
            </>
          )}
        </p>
      ) : (
        <p className="rg-job-card__note">The pull was cancelled.</p>
      )}

      {!active && (
        <footer className="rg-job-card__actions">
          <button
            type="button"
            className="rg-btn rg-btn--sm rg-btn--bare"
            onClick={() => dismiss(entry.id)}
          >
            Dismiss
          </button>
        </footer>
      )}
    </article>
  );
}

export function BridgeJobsDock() {
  useBridgeJobPolling();
  const bridge = useBridge();
  const entries = useBridgeJobEntries();
  const clearFinished = useBridgeJobStore((state) => state.clearFinished);
  const [collapsed, setCollapsed] = useState(false);
  const known = useRef<Set<string>>(new Set());

  // A newly started pull re-opens a collapsed dock, as the local console does.
  useEffect(() => {
    let appeared = false;
    for (const entry of entries) {
      if (!known.current.has(entry.id)) {
        known.current.add(entry.id);
        appeared = true;
      }
    }
    if (appeared) setCollapsed(false);
  }, [entries]);

  if (entries.length === 0) return null;

  const activeCount = entries.filter(isActiveEntry).length;
  const finishedCount = entries.length - activeCount;
  const instanceId = bridge.ping?.instance_id ?? '';

  return (
    <aside
      className={cn('rg-job-console', collapsed && 'rg-job-console--collapsed')}
      aria-label="Local pulls"
    >
      <header className="rg-job-console__header">
        <button
          type="button"
          className="rg-job-console__toggle"
          onClick={() => setCollapsed((current) => !current)}
          aria-expanded={!collapsed}
        >
          <span aria-hidden="true">{collapsed ? '▲' : '▼'}</span> Pulls to your local refgenie
        </button>
        <span className="rg-job-console__badge" aria-live="polite">
          {activeCount} active
        </span>
        <span className="rg-job-console__status rg-job-console__status--open">
          <span className="rg-job-console__dot" aria-hidden="true" />
          {bridge.ping?.instance_label}
        </span>
        <Link className="rg-link text-sm" to="/local/jobs">
          All local jobs
        </Link>
        {finishedCount > 0 && (
          <button
            type="button"
            className="rg-btn rg-btn--sm rg-btn--bare"
            onClick={() => clearFinished(instanceId)}
          >
            Clear finished
          </button>
        )}
      </header>

      {!collapsed && (
        <div className="rg-job-console__body">
          {entries.map((entry) => (
            <BridgeJobCard key={entry.id} entry={entry} localBaseUrl={bridge.baseUrl} />
          ))}
        </div>
      )}
    </aside>
  );
}
