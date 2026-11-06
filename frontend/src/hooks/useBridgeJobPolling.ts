/**
 * Polls the connected local refgenie for each active bridge-started pull,
 * mounted once by `BridgeJobsDock`. It stops when the bridge disconnects, and
 * reports completion here so the toast still fires after the user navigates away.
 */

import { useEffect, useMemo } from 'react';
import { getJob } from '../services/jobs';
import { isTerminal } from '../services/contracts';
import { ApiError } from '../services/http';
import { createBridgeClients } from '../services/bridge/client';
import { BRIDGE_POLL_INTERVAL_MS } from '../services/bridge/pull';
import { useBridgeStore } from '../stores/bridgeStore';
import { entriesFor, isActiveEntry, useBridgeJobStore } from '../stores/bridgeJobStore';
import type { BridgeJobEntry } from '../stores/bridgeJobStore';
import { useToast } from './useToast';
import { useInvalidateBridgeDigests } from './queries/useBridgeDigests';

export function useBridgeJobPolling(intervalMs = BRIDGE_POLL_INTERVAL_MS): void {
  const status = useBridgeStore((s) => s.status);
  const baseUrl = useBridgeStore((s) => s.baseUrl);
  const instanceId = useBridgeStore((s) => s.ping?.instance_id ?? null);
  // A string, so the effect restarts only when the SET of active ids changes,
  // not on every progress update.
  const activeIds = useBridgeJobStore((s) =>
    entriesFor(s.entries, instanceId)
      .filter(isActiveEntry)
      .map((entry) => entry.id)
      .join(','),
  );
  const toast = useToast();
  const invalidate = useInvalidateBridgeDigests();

  useEffect(() => {
    if (status !== 'connected' || !baseUrl || !activeIds) return;
    const controller = new AbortController();
    const { local } = createBridgeClients(baseUrl);
    const ids = activeIds.split(',');
    let timer: ReturnType<typeof setTimeout> | undefined;

    const pollOne = async (id: string) => {
      try {
        const job = await getJob(local, id, { signal: controller.signal });
        if (controller.signal.aborted) return;
        const store = useBridgeJobStore.getState();
        const before = store.entries.find((entry) => entry.id === id);
        store.update(id, job);
        if (!before || !isActiveEntry(before) || !isTerminal(job.status)) return;
        const label = before.label;
        if (job.status === 'succeeded') {
          toast.success(`${label} is now on your local refgenie.`);
          void invalidate();
        } else if (job.status === 'failed') {
          toast.error(job.error?.message ?? `Pulling ${label} failed.`);
        }
      } catch (error) {
        if (controller.signal.aborted) return;
        // Jobs are not persisted by the local server: after a restart it no
        // longer knows the id. Anything else is transient; try again next tick.
        if (error instanceof ApiError && error.status === 404) {
          useBridgeJobStore.getState().markLost(id);
        }
      }
    };

    const tick = async () => {
      await Promise.all(ids.map(pollOne));
      if (!controller.signal.aborted) timer = setTimeout(() => void tick(), intervalMs);
    };
    void tick();

    return () => {
      controller.abort();
      if (timer) clearTimeout(timer);
    };
  }, [status, baseUrl, activeIds, intervalMs, toast, invalidate]);
}

/** Tracked pulls of the currently connected local refgenie; empty while disconnected. */
export function useBridgeJobEntries(): BridgeJobEntry[] {
  const status = useBridgeStore((s) => s.status);
  const instanceId = useBridgeStore((s) => s.ping?.instance_id ?? null);
  const entries = useBridgeJobStore((s) => s.entries);
  return useMemo(
    () => (status === 'connected' ? entriesFor(entries, instanceId) : []),
    [status, entries, instanceId],
  );
}
