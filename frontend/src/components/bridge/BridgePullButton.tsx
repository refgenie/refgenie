/**
 * "Pull to my refgenie" (bridge plan B6) — the CROSS-ORIGIN pull, distinct
 * from `components/actions/PullButton.tsx`, which pulls on this instance.
 *
 * Rendered only on the remote branch with the bridge connected. When
 * cross-origin actions are disabled on the local side
 * (`bridge.actions_cross_origin` false — the default, bridge mode `read`), the
 * button is disabled with the remedy in its title and a deep link into the
 * local SPA's prefilled, UNSUBMITTED pull page is offered instead. That page
 * never auto-executes, here or anywhere. The same deep link is offered when
 * the local side refuses the pull with `server_not_subscribed`: a bridge page
 * may only pull from servers the user subscribed to, so any other server has
 * to be confirmed on the local origin.
 *
 * The job runs on a machine `jobStore` does not manage, so it is deliberately
 * not registered there and `usePullAction` is not reused. It is registered in
 * `bridgeJobStore` instead: the app-level poller and `BridgeJobsDock` report
 * it from then on, so progress survives navigation and reload, and this
 * button only mirrors the active job for its own target.
 */

import { useState } from 'react';
import { ExternalLink } from '../common/ExternalLink';
import { JobProgress } from '../jobs/JobProgress';
import { phaseInfo } from '../jobs/phases';
import { useApiClient } from '../../hooks/useApiClient';
import { useBridge } from '../../hooks/useBridge';
import { useRouteBase } from '../../hooks/useRouteBase';
import { useToast } from '../../hooks/useToast';
import { absoluteServerRoot } from '../../services/instanceLinks';
import { createBridgeClients } from '../../services/bridge/client';
import { pullAsset } from '../../services/actions';
import { ApiError } from '../../services/http';
import { isActiveEntry, useBridgeJobStore } from '../../stores/bridgeJobStore';

export interface BridgePullButtonProps {
  genomeDigest: string;
  assetGroupName: string;
  assetName?: string | null;
  /**
   * Shorten the labels for a table cell, where the full sentence repeats on
   * every row and the surrounding column already says what the row is.
   */
  compact?: boolean;
}

export function BridgePullButton({
  genomeDigest,
  assetGroupName,
  assetName,
  compact = false,
}: BridgePullButtonProps) {
  const bridge = useBridge();
  const client = useApiClient(); // the REMOTE read client
  const routeBase = useRouteBase();
  const toast = useToast();
  const [pending, setPending] = useState(false);
  // Set when the local refgenie refused because it does not subscribe to this
  // server: the user has to confirm the pull on their own origin instead.
  const [needsLocalConfirm, setNeedsLocalConfirm] = useState(false);
  const instanceId = bridge.ping?.instance_id ?? null;
  const activeJob = useBridgeJobStore((state) =>
    state.entries.find(
      (entry) =>
        entry.instanceId === instanceId &&
        entry.genomeDigest === genomeDigest &&
        entry.assetGroupName === assetGroupName &&
        entry.assetName === (assetName ?? null) &&
        isActiveEntry(entry),
    ),
  );

  // Rendered ONLY on the remote branch with the bridge connected, and gated on
  // the LOCAL instance's pull capability from /ping. `useCapability('pull')`
  // would be wrong: on refgenie.org the page's own capability is false, which
  // would silently kill the whole feature.
  if (
    bridge.status !== 'connected' ||
    !bridge.baseUrl ||
    routeBase !== '' ||
    bridge.ping?.capabilities.pull !== true
  ) {
    return null;
  }

  const baseUrl = bridge.baseUrl;
  // `server_url` is the remote server ROOT, and it MUST be absolute: it is read
  // by the local refgenie, on its own origin, where the same-origin `/v4`
  // deployment's bare `serverRoot` of `''` addresses nothing. Sending that
  // empty string silently drops the "pull from THIS server" instruction and
  // leaves the local instance falling back to its own subscriptions.
  const remoteBase = absoluteServerRoot(client.baseUrl);
  const crossOriginEnabled = bridge.ping?.bridge.actions_cross_origin === true;
  const registryLabel = `${assetGroupName}${assetName ? `/${assetName}` : ''}`;
  const deepLink =
    `${baseUrl}/pull?server=${encodeURIComponent(remoteBase)}` +
    `&genome=${encodeURIComponent(genomeDigest)}` +
    `&asset_group=${encodeURIComponent(assetGroupName)}` +
    (assetName ? `&asset=${encodeURIComponent(assetName)}` : '');

  const startPull = async () => {
    setPending(true);
    try {
      const { local } = createBridgeClients(baseUrl);
      const ref = await pullAsset(local, {
        asset_group: assetGroupName,
        // Exactly one genome reference: `genome` is omitted entirely, not
        // nulled, because both-or-neither is a 422.
        genome_digest: genomeDigest,
        asset: assetName ?? null,
        server_url: remoteBase,
        // ALWAYS explicit — see frontend/README.md. Omitting it lets the
        // puller reach a stdin prompt that hangs a worker forever.
        force: false,
      });
      useBridgeJobStore.getState().track({
        id: ref.job_id,
        instanceId: bridge.ping?.instance_id ?? '',
        label: registryLabel,
        genomeDigest,
        assetGroupName,
        assetName: assetName ?? null,
        submittedAt: new Date().toISOString(),
      });
    } catch (error) {
      if (error instanceof ApiError && error.code === 'server_not_subscribed') {
        setNeedsLocalConfirm(true);
        toast.error(error.detail);
        return;
      }
      toast.error(
        error instanceof ApiError ? error.detail : 'Could not reach the local refgenie.',
      );
    } finally {
      setPending(false);
    }
  };

  if (activeJob) {
    const progress = activeJob.job?.status === 'running' ? activeJob.job.progress : null;
    return (
      <span className="rg-pull-button rg-pull-button--active">
        <JobProgress kind="pull" progress={progress} compact />
        <span className="rg-pull-button__phase rg-muted">
          {activeJob.job?.status === 'running'
            ? phaseInfo('pull', progress?.phase).label
            : 'Queued'}
        </span>
      </span>
    );
  }

  return (
    <span className="rg-pull-button">
      <button
        type="button"
        className="rg-btn rg-btn--sm rg-btn--primary"
        disabled={!crossOriginEnabled || pending}
        title={
          crossOriginEnabled
            ? compact
              ? `Pull ${registryLabel} to your local refgenie`
              : undefined
            : `Cross-origin actions are disabled on this refgenie (bridge mode is '${bridge.ping?.bridge_mode}'). Restart with \`refgenie dash --bridge full\` to enable them.`
        }
        onClick={() => void startPull()}
      >
        {pending ? 'Pulling…' : compact ? 'Pull' : 'Pull to my refgenie'}
      </button>
      {(!crossOriginEnabled || needsLocalConfirm) && (
        <ExternalLink
          className="rg-btn rg-btn--sm"
          href={deepLink}
          title="Opens your local refgenie with this pull prefilled (nothing runs until you confirm there)."
        >
          {compact ? 'Open locally' : 'Open in local refgenie'}
        </ExternalLink>
      )}
    </span>
  );
}
