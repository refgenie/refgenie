/**
 * Silently re-probes a remembered local connection on later visits. The first
 * connect is always user-initiated. Mounted once by AppLayout.
 *
 * The re-probe runs only when it cannot raise the browser's local-network
 * permission prompt: a visitor who connected once should not be asked on
 * every visit. When the browser would prompt, the user reconnects by clicking.
 */

import { useEffect } from 'react';
import { useBridge } from './useBridge';
import { isSelfLocal, localNetworkPermission } from '../services/bridge/environment';
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
      void localNetworkPermission().then((permission) => {
        if (permission !== 'granted' && permission !== 'unsupported') return;
        if (useBridgeStore.getState().status !== 'idle') return;
        void connect(remembered.port, { silent: true });
      });
    } else if (isSelfLocal(state.port)) {
      useBridgeStore.getState().setStatus('self');
    }
  }, [enabled, connect]);
}
