/**
 * Bridge state plus connect and disconnect. Has no effects, since many components
 * mount it; the auto re-probe lives in `useBridgeAutoConnect`.
 */

import { useCallback } from 'react';
import { isSelfLocal, isWebKitOnly } from '../services/bridge/environment';
import { forget, remember } from '../services/bridge/persistence';
import { probeLocal } from '../services/bridge/probe';
import { useBridgeStore } from '../stores/bridgeStore';
import { useBridgeJobStore } from '../stores/bridgeJobStore';

export interface BridgeConnectOptions {
  /** A silent re-probe never renders the troubleshooting panel. */
  silent?: boolean;
}

export function useBridge() {
  const store = useBridgeStore();

  const connect = useCallback(async (port: number, options?: BridgeConnectOptions) => {
    const silent = options?.silent ?? false;
    const state = useBridgeStore.getState();
    if (isSelfLocal(port)) {
      // This page is the local dash itself.
      state.setStatus('self');
      return;
    }
    if (isWebKitOnly()) {
      // WebKit blocks the request outright, so skip it.
      state.setStatus('blocked');
      return;
    }
    state.setPort(port);
    state.setConnecting();
    const result = await probeLocal(port);
    const current = useBridgeStore.getState();
    if (result.kind === 'connected') {
      current.setConnected(result.ping, port, result.baseUrl);
      remember({
        port,
        instanceId: result.ping.instance_id,
        connectedAt: new Date().toISOString(),
      });
    } else if (result.kind === 'unsupported') {
      current.setUnsupported(result.reason);
    } else {
      // Show troubleshooting only after an explicit click.
      current.setAbsent(!silent);
    }
  }, []);

  const disconnect = useCallback(() => {
    forget();
    // Forget the pulls tracked on this connection too.
    useBridgeJobStore.getState().clear();
    useBridgeStore.getState().reset();
  }, []);

  return { ...store, connect, disconnect };
}
