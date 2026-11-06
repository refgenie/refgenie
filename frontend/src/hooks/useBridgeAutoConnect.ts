/**
 * Silently re-probes a remembered local connection on later visits. The first
 * connect is always user-initiated. Mounted once by AppLayout.
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
