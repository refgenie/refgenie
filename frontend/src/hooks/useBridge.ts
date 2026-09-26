/**
 * The bridge's store read plus its two commands.
 *
 * Deliberately effect-free: every component that shows bridge state mounts
 * this, and a probe fired from here would fire once per consumer. The one
 * remembered-connection re-probe lives in `useBridgeAutoConnect`, mounted once
 * by AppLayout.
 */

import { useCallback } from 'react';
import { isSelfLocal, isWebKitOnly } from '../services/bridge/environment';
import { forget, remember } from '../services/bridge/persistence';
import { probeLocal } from '../services/bridge/probe';
import { useBridgeStore } from '../stores/bridgeStore';

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
      // D7 second gate; the first is `config.mode !== 'local'`, checked before
      // any bridge component mounts.
      state.setStatus('self');
      return;
    }
    if (isWebKitOnly()) {
      // Skip the doomed request: WebKit blocks it outright, and this is the
      // one failure cause identifiable with confidence.
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
      // Never nag on a normal page load: the troubleshooting panel renders
      // only after an explicit click.
      current.setAbsent(!silent);
    }
  }, []);

  const disconnect = useCallback(() => {
    forget();
    useBridgeStore.getState().reset();
  }, []);

  return { ...store, connect, disconnect };
}
