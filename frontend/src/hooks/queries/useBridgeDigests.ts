/**
 * The connected local refgenie's digest index, and its invalidator.
 *
 * The read is scoped to whichever cache the tree is under, which is always the
 * root one in practice: presence badges render on the REMOTE branch, against
 * the page's own listing, and callers inside `/local` pass `enabled: false`.
 */

import { useMemo } from 'react';
import { useInvalidateResources, useResource } from '../useResource';
import { qk } from '../../services/queryKeys';
import { createBridgeClients } from '../../services/bridge/client';
import { fetchBridgeDigestIndex } from '../../services/bridge/digestIndex';
import { useBridgeStore } from '../../stores/bridgeStore';

/**
 * @param enabled Callers on the `/local` branch pass `false`: inside that scope
 *   every row is local by definition, so the index would be a wasted request.
 */
export function useBridgeDigests(enabled = true) {
  const status = useBridgeStore((s) => s.status);
  const baseUrl = useBridgeStore((s) => s.baseUrl);
  const clients = useMemo(
    () => (baseUrl ? createBridgeClients(baseUrl) : null),
    [baseUrl],
  );
  return useResource(
    qk.bridgeDigests(baseUrl ?? ''),
    () => fetchBridgeDigestIndex(clients!.read),
    {
      enabled: enabled && status === 'connected' && !!clients,
      staleTime: 60_000,
    },
  );
}

/**
 * Invalidate after a successful bridge pull so badges update immediately.
 *
 * Deliberately NOT published as a domain key on `services/invalidation.ts`:
 * that bus invalidates *this* backend's reads, and a pull that landed on a
 * different machine must not touch them.
 */
export function useInvalidateBridgeDigests() {
  const invalidate = useInvalidateResources();
  return useMemo(() => () => invalidate(['bridge-digests']), [invalidate]);
}
