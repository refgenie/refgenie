/**
 * Waiting on a job that is running on a DIFFERENT refgenie.
 *
 * Submission itself is `services/actions.ts::pullAsset` — the bridge target is
 * a refgenie1 server, so it shares the request model, the paths, the error
 * envelope and the `X-Refgenie-Action` header with the local surface. The only
 * thing that cannot be shared is progress: `useJobEvents` streams THIS
 * instance's jobs into `jobStore`, and a cross-origin instance's jobs must
 * never enter that store or the console and the invalidation bus would both be
 * reporting on a machine they do not manage. So the bridge polls, from ONE
 * place (`hooks/useBridgeJobPolling.ts`) over the jobs remembered in
 * `stores/bridgeJobStore.ts`, so progress outlives the button that started it.
 *
 * A cross-origin SSE stream to `/v1/jobs/events` would technically work
 * (EventSource sends no headers, the endpoint is an unguarded GET, and bridge
 * origins are allowed). It is deliberately not used: a long-lived connection
 * from a public page into the user's machine is a larger commitment than a
 * poll that terminates.
 */

/** How often `hooks/useBridgeJobPolling.ts` re-reads each active bridge job. */
export const BRIDGE_POLL_INTERVAL_MS = 1000;
