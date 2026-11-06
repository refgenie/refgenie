/** Localhost-bridge connection state, shared by components in unrelated subtrees. */

import { create } from 'zustand';
import { DEFAULT_LOCAL_PORT } from '../services/bridge/persistence';
import type { PingResponse } from '../services/bridge/contract';

/**
 * - `idle`: not probed yet
 * - `connecting`: probe in flight
 * - `connected`: a supported local refgenie answered
 * - `absent`: probe failed (not running, access denied, wrong port; a page cannot tell which)
 * - `blocked`: WebKit-only browser; never probed
 * - `unsupported`: something answered, but not a bridge we support
 * - `self`: this page is the local dash; hide connect UI
 */
export type BridgeStatus =
  | 'idle'
  | 'connecting'
  | 'connected'
  | 'absent'
  | 'blocked'
  | 'unsupported'
  | 'self';

export type UnsupportedReason =
  | 'wrong-service'
  | 'unsupported-version'
  | null;

export interface BridgeStoreState {
  status: BridgeStatus;
  /** The port the next probe targets. Survives `reset` so a reconnect keeps it. */
  port: number;
  /** `http://localhost:{port}` while connected, else null. */
  baseUrl: string | null;
  ping: PingResponse | null;
  unsupportedReason: UnsupportedReason;
  /** True when the last failed probe was a user click, so troubleshooting help shows. */
  showTroubleshooting: boolean;
  /** One auto-probe per page load. Kept in the store, not a ref, so StrictMode cannot probe twice. */
  autoProbeDone: boolean;

  setPort: (port: number) => void;
  setStatus: (status: BridgeStatus) => void;
  setConnecting: () => void;
  setConnected: (ping: PingResponse, port: number, baseUrl: string) => void;
  setAbsent: (explicit: boolean) => void;
  setUnsupported: (reason: Exclude<UnsupportedReason, null>) => void;
  markAutoProbeDone: () => void;
  reset: () => void;
}

/** Connection state that `reset` clears. Excludes `port` (kept across reconnects) and `autoProbeDone`. */
const EMPTY: Pick<
  BridgeStoreState,
  'status' | 'baseUrl' | 'ping' | 'unsupportedReason' | 'showTroubleshooting'
> = {
  status: 'idle',
  baseUrl: null,
  ping: null,
  unsupportedReason: null,
  showTroubleshooting: false,
};

export const useBridgeStore = create<BridgeStoreState>()((set) => ({
  ...EMPTY,
  port: DEFAULT_LOCAL_PORT,
  autoProbeDone: false,

  setPort: (port) => set({ port }),

  setStatus: (status) => set({ status }),

  setConnecting: () => set({ status: 'connecting', showTroubleshooting: false }),

  setConnected: (ping, port, baseUrl) =>
    set({
      status: 'connected',
      ping,
      port,
      baseUrl,
      unsupportedReason: null,
      showTroubleshooting: false,
    }),

  setAbsent: (explicit) =>
    set({
      status: 'absent',
      baseUrl: null,
      ping: null,
      showTroubleshooting: explicit,
    }),

  setUnsupported: (reason) =>
    set({
      status: 'unsupported',
      baseUrl: null,
      ping: null,
      unsupportedReason: reason,
    }),

  markAutoProbeDone: () => set({ autoProbeDone: true }),

  reset: () => set({ ...EMPTY }),
}));
