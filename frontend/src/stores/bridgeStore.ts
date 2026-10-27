/**
 * Localhost-bridge connection state.
 *
 * zustand rather than React context for the same reason `jobStore` is: the
 * bridge is read by the sidebar control, the connect card, the connect dialog,
 * the presence badges and the pull affordance, none of which share a subtree.
 *
 * The store is a plain module with no React import beyond the hook factory, so
 * `useBridge` can drive it and tests can exercise it without rendering.
 */

import { create } from 'zustand';
import { DEFAULT_LOCAL_PORT } from '../services/bridge/persistence';
import type { PingResponse } from '../services/bridge/contract';

/**
 * - `idle`        — never probed this visit, no remembered connection tried
 * - `connecting`  — a probe is in flight
 * - `connected`   — a supported local refgenie answered
 * - `absent`      — the probe failed (any of: not running, LNA denied,
 *                   wrong port, bridge off — indistinguishable from a page)
 * - `blocked`     — WebKit-only browser; the probe is never attempted
 * - `unsupported` — something answered but is not a bridge we can speak to
 * - `self`        — this page IS the local dash (D7); hide all connect UI
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
  /**
   * True when the last failed probe was an explicit user click (render the
   * troubleshooting panel); false for silent remembered re-probes.
   */
  showTroubleshooting: boolean;
  /**
   * One auto-probe attempt per page load, tracked here rather than in a ref so
   * StrictMode's double effect invocation cannot fire two probes at the user's
   * machine.
   */
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

/**
 * The connection-shaped half of the state. `port` is excluded deliberately (a
 * disconnect keeps the port the user typed), and so is `autoProbeDone`, which
 * is a per-page-load latch rather than connection state.
 */
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
