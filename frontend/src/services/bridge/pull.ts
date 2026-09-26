/**
 * Waiting on a job that is running on a DIFFERENT refgenie.
 *
 * Submission itself is `services/actions.ts::pullAsset` — the bridge target is
 * a refgenie1 server, so it shares the request model, the paths, the error
 * envelope and the `X-Refgenie-Action` header with the local surface. The only
 * thing that cannot be shared is progress: `useJobEvents` streams THIS
 * instance's jobs into `jobStore`, and a cross-origin instance's jobs must
 * never enter that store or the console and the invalidation bus would both be
 * reporting on a machine they do not manage. So the bridge polls.
 *
 * A cross-origin SSE stream to `/v1/jobs/events` would technically work
 * (EventSource sends no headers, the endpoint is an unguarded GET, and bridge
 * origins are allowed). It is deliberately not used: a long-lived connection
 * from a public page into the user's machine is a larger commitment than a
 * poll that terminates.
 */

import { getJob } from '../jobs';
import { isTerminal } from '../contracts';
import type { Job } from '../contracts';
import type { ApiClient } from '../http';

export const BRIDGE_POLL_INTERVAL_MS = 1000;

/**
 * Poll until the job is terminal; `onUpdate` fires on every poll. The signal
 * is not optional in practice: without it an unmounted button keeps polling
 * the user's machine for the life of the tab.
 */
export async function waitForBridgeJob(
  client: ApiClient,
  jobId: string,
  onUpdate: (job: Job) => void,
  signal: AbortSignal,
  intervalMs = BRIDGE_POLL_INTERVAL_MS,
): Promise<Job | null> {
  for (;;) {
    if (signal.aborted) return null;
    const job = await getJob(client, jobId, { signal });
    onUpdate(job);
    if (isTerminal(job.status)) return job;
    await new Promise((resolve) => setTimeout(resolve, intervalMs));
  }
}
