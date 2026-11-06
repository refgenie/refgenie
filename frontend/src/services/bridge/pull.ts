/**
 * Bridge pulls run on another refgenie, so their progress is polled (in
 * `hooks/useBridgeJobPolling.ts`) rather than streamed into this instance's
 * `jobStore`. Polling also avoids holding a long-lived connection into the
 * user's machine.
 */

/** How often `hooks/useBridgeJobPolling.ts` re-reads each active bridge job. */
export const BRIDGE_POLL_INTERVAL_MS = 1000;
