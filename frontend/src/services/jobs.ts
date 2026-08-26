/**
 * The jobs service: reads and the one SSE URL.
 *
 * Two constraints are baked into the shape here and must not be relaxed:
 *  - **One stream, not one per job.** `jobEventsUrl` takes no job id. Browsers
 *    cap ~6 HTTP/1.1 connections per origin, and a user watching four pulls
 *    would exhaust that.
 *  - **No custom header on the stream.** `EventSource` cannot set request
 *    headers, so the events endpoint is a plain GET; `X-Refgenie-Action`
 *    applies only to state-changing calls, which go through `mutate`.
 */

import { JOB_PATHS } from './contracts';
import type { Job, JobLogPage, JobStatus } from './contracts';
import type { ApiClient, RequestInitLite } from './http';
import type { Paginated } from '../types/pagination';

/** `active` is queued+running; `terminal` is succeeded/failed/cancelled. */
export type JobStatusFilter = 'active' | 'terminal' | JobStatus;

export interface ListJobsParams {
  status?: JobStatusFilter;
  kind?: string;
  offset?: number;
  limit?: number;
}

/** The seam `useJobEvents` polls through when SSE is unavailable. */
export type JobsFetcher = (
  params: ListJobsParams,
  init?: RequestInitLite,
) => Promise<Paginated<Job>>;

export const listJobs = (
  c: ApiClient,
  p: ListJobsParams = {},
  init?: RequestInitLite,
): Promise<Paginated<Job>> =>
  c.get<Paginated<Job>>(
    JOB_PATHS.list,
    { status: p.status, kind: p.kind, offset: p.offset, limit: p.limit },
    init,
  );

export const getJob = (c: ApiClient, id: string, init?: RequestInitLite): Promise<Job> =>
  c.get<Job>(JOB_PATHS.job(id), undefined, init);

export const getJobLog = (
  c: ApiClient,
  id: string,
  offset = 0,
  init?: RequestInitLite,
): Promise<JobLogPage> => c.get<JobLogPage>(JOB_PATHS.log(id), { offset }, init);

/** Gated on `capabilities.jobs_cancel`; a running build answers 409. */
export const cancelJob = (c: ApiClient, id: string, init?: RequestInitLite): Promise<void> =>
  c.mutate<void>(JOB_PATHS.cancel(id), {
    method: 'POST',
    action: 'job.cancel',
    signal: init?.signal,
  });

/**
 * The multiplexed stream URL. `since` resumes from a manager-global seq, so a
 * reconnect replays whatever happened while disconnected instead of losing it.
 */
export function jobEventsUrl(c: ApiClient, since?: number): string {
  return c.url(JOB_PATHS.events, since && since > 0 ? { since } : undefined);
}

/** Build a fetcher bound to a client, for the hook's polling fallback. */
export function jobsFetcher(c: ApiClient): JobsFetcher {
  return (params, init) => listJobs(c, params, init);
}
