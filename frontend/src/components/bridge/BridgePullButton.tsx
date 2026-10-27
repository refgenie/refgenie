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
 * not registered there and `usePullAction` is not reused: this polls the other
 * instance and reports inline.
 */

import { useEffect, useRef, useState } from 'react';
import { ExternalLink } from '../common/ExternalLink';
import { JobProgress } from '../jobs/JobProgress';
import { phaseInfo } from '../jobs/phases';
import { useApiClient } from '../../hooks/useApiClient';
import { useBridge } from '../../hooks/useBridge';
import { useRouteBase } from '../../hooks/useRouteBase';
import { useToast } from '../../hooks/useToast';
import { useInvalidateBridgeDigests } from '../../hooks/queries/useBridgeDigests';
import { absoluteServerRoot } from '../../services/instanceLinks';
import { createBridgeClients } from '../../services/bridge/client';
import { waitForBridgeJob } from '../../services/bridge/pull';
import { pullAsset } from '../../services/actions';
import { ApiError } from '../../services/http';
import { isTerminal } from '../../services/contracts';
import type { Job } from '../../services/contracts';

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
  const invalidate = useInvalidateBridgeDigests();
  const [job, setJob] = useState<Job | null>(null);
  const [pending, setPending] = useState(false);
  // Set when the local refgenie refused because it does not subscribe to this
  // server: the user has to confirm the pull on their own origin instead.
  const [needsLocalConfirm, setNeedsLocalConfirm] = useState(false);
  const aborter = useRef<AbortController | null>(null);

  // Unmounting must stop the poll; without this the button keeps hitting the
  // user's machine for the life of the tab.
  useEffect(() => () => aborter.current?.abort(), []);

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
    const controller = new AbortController();
    aborter.current?.abort();
    aborter.current = controller;
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
      const final = await waitForBridgeJob(
        local,
        ref.job_id,
        setJob,
        controller.signal,
      );
      if (!final) return;
      if (final.status === 'succeeded') {
        toast.success(`${registryLabel} is now on your local refgenie.`);
        void invalidate();
      } else {
        toast.error(
          final.error?.message ?? `Job ended with status '${final.status}'.`,
        );
      }
    } catch (error) {
      if (controller.signal.aborted) return;
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
      setJob(null);
    }
  };

  if (job && !isTerminal(job.status)) {
    return (
      <span className="rg-pull-button rg-pull-button--active">
        <JobProgress kind="pull" progress={job.progress} compact />
        <span className="rg-pull-button__phase rg-muted">
          {phaseInfo('pull', job.progress?.phase).label}
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
