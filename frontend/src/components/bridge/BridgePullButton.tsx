/**
 * Pulls a remote asset into the user's local refgenie over the bridge. If the
 * local side disallows cross-origin actions or is not subscribed to this
 * server, it links to the local SPA's prefilled pull page instead. Jobs are
 * tracked in `bridgeJobStore`, not `jobStore`, because they run on another machine.
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
import { localServersUrl } from '../servers/links';
import { ApiError } from '../../services/http';
import { isActiveEntry, useBridgeJobStore } from '../../stores/bridgeJobStore';

export interface BridgePullButtonProps {
  genomeDigest: string;
  assetGroupName: string;
  assetName?: string | null;
  /** Use short labels, for table cells. */
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
  // Set when the local refgenie is not subscribed to this server.
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

  // Gate on the LOCAL instance's pull capability from /ping, not
  // `useCapability('pull')`, which reflects this (remote) page's server.
  if (
    bridge.status !== 'connected' ||
    !bridge.baseUrl ||
    routeBase !== '' ||
    bridge.ping?.capabilities.pull !== true
  ) {
    return null;
  }

  const baseUrl = bridge.baseUrl;
  // Must be absolute: the local refgenie reads it on its own origin, where an
  // empty same-origin root would make it fall back to its own subscriptions.
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
        // Send exactly one genome reference; `genome` is omitted, not nulled (422 otherwise).
        genome_digest: genomeDigest,
        asset: assetName ?? null,
        server_url: remoteBase,
        // Always explicit: without it the puller can hang on a stdin prompt.
        force: false,
      });
      useBridgeJobStore.getState().track({
        id: ref.job_id,
        instanceId: bridge.ping?.instance_id ?? '',
        label: registryLabel,
        genomeDigest,
        assetGroupName,
        assetName: assetName ?? null,
        serverUrl: remoteBase,
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
      {needsLocalConfirm && (
        // Subscribing is same-origin only, so link to the local prefilled form.
        <ExternalLink
          className="rg-btn rg-btn--sm"
          href={localServersUrl(baseUrl, remoteBase)}
          title="Opens your local refgenie's Servers page with this server filled in (nothing changes until you click Subscribe there)."
        >
          {compact ? 'Subscribe locally' : 'Subscribe to this server locally'}
        </ExternalLink>
      )}
    </span>
  );
}
