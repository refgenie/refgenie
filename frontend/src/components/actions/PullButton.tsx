/**
 * Pull button on a remote asset row. Shows the active job's progress while a
 * matching pull is queued or running.
 */

import { useState } from 'react';
import { useCapability } from '../../hooks/useCapability';
import { usePullAction } from '../../hooks/usePullAction';
import { targetKey, useJobStore } from '../../stores/jobStore';
import { JobProgress } from '../jobs/JobProgress';
import { phaseInfo } from '../jobs/phases';

export interface PullButtonProps {
  serverUrl: string | null;
  genomeDigest: string;
  genomeName?: string | null;
  assetGroupName: string;
  assetName: string;
  existsLocally: boolean;
  /** Offered by the "try another server" path; omitted when there is one server. */
  onPickServer?: () => void;
}

export function PullButton({
  serverUrl,
  genomeDigest,
  genomeName,
  assetGroupName,
  assetName,
  existsLocally,
  onPickServer,
}: PullButtonProps) {
  const canPull = useCapability('pull');
  const { pull, isPending } = usePullAction();
  const [menuOpen, setMenuOpen] = useState(false);

  const key = targetKey('pull', {
    genome_digest: genomeDigest,
    genome_name: genomeName ?? null,
    asset_group_name: assetGroupName,
    asset_name: assetName,
  });
  const active = useJobStore((state) => state.activeJobFor(key));
  const focusJob = useJobStore((state) => state.focusJob);

  if (!canPull) return null;

  const submit = (force: boolean) => {
    setMenuOpen(false);
    void pull({
      genomeDigest,
      genomeName,
      assetGroupName,
      assetName,
      serverUrl,
      force,
    });
  };

  if (active) {
    return (
      <span className="rg-pull-button rg-pull-button--active">
        <JobProgress kind={active.kind} progress={active.progress} compact />
        <span className="rg-pull-button__phase rg-muted">
          {phaseInfo(active.kind, active.progress?.phase).label}
        </span>
        <button
          type="button"
          className="rg-btn rg-btn--sm rg-btn--bare"
          onClick={() => focusJob(active.id)}
        >
          Show console
        </button>
      </span>
    );
  }

  return (
    <span className="rg-pull-button">
      <button
        type="button"
        className="rg-btn rg-btn--sm rg-btn--primary"
        disabled={existsLocally || isPending}
        onClick={() => submit(false)}
      >
        {existsLocally ? 'Already local' : 'Pull'}
      </button>

      <button
        type="button"
        className="rg-btn rg-btn--sm"
        aria-label="More pull options"
        aria-expanded={menuOpen}
        onClick={() => setMenuOpen((open) => !open)}
      >
        &hellip;
      </button>

      {menuOpen && (
        <span className="rg-pull-button__menu" role="menu">
          <button
            type="button"
            className="rg-btn rg-btn--sm"
            role="menuitem"
            onClick={() => submit(true)}
          >
            {existsLocally ? 'Re-pull (overwrite)' : 'Pull anyway (overwrite)'}
          </button>
          {onPickServer && (
            <button
              type="button"
              className="rg-btn rg-btn--sm"
              role="menuitem"
              onClick={() => {
                setMenuOpen(false);
                onPickServer();
              }}
            >
              Try another server
            </button>
          )}
        </span>
      )}
    </span>
  );
}
