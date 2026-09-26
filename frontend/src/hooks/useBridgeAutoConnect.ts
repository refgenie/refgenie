/**
 * The one remembered-connection re-probe, mounted once by AppLayout alongside
 * useJobEvents and useInvalidationBridge.
 *
 * Detection is user-initiated on first connect and remembered afterward (D3):
 * no silent probe on a cold visit — that is both the Chrome-LNA-friendly
 * behavior and the honest one. Once a connection has succeeded, the remembered
 * `{port, instanceId}` licenses an automatic, silent re-probe on later visits.
 */

import { useEffect } from 'react';
import { useBridge } from './useBridge';
import { isSelfLocal } from '../services/bridge/environment';
import { recall } from '../services/bridge/persistence';
import { useBridgeStore } from '../stores/bridgeStore';

export function useBridgeAutoConnect(enabled: boolean) {
  const { connect } = useBridge();

  useEffect(() => {
    if (!enabled) return;
    const state = useBridgeStore.getState();
    if (state.autoProbeDone) return;
    state.markAutoProbeDone();
    const remembered = recall();
    if (remembered && state.status === 'idle') {
      void connect(remembered.port, { silent: true });
    } else if (isSelfLocal(state.port)) {
      useBridgeStore.getState().setStatus('self');
    }
  }, [enabled, connect]);
}
