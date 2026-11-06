/** Digest index of the connected local refgenie, used for presence badges, and its invalidator. */

import { useMemo } from 'react';
import { useInvalidate, useResource } from '../useResource';
import { qk } from '../../services/queryKeys';
import { createBridgeClients } from '../../services/bridge/client';
import { fetchBridgeDigestIndex } from '../../services/bridge/digestIndex';
import { useBridgeStore } from '../../stores/bridgeStore';

/** @param enabled Pass `false` under `/local`, where every row is local already. */
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

/** Invalidate after a successful bridge pull so badges update immediately. */
export function useInvalidateBridgeDigests() {
  const invalidate = useInvalidate();
  return useMemo(() => () => invalidate(['bridge-digests']), [invalidate]);
}
