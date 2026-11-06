/**
 * Job reads and the event stream URL. There is one stream for all jobs, because
 * browsers allow only ~6 connections per origin.
 */

import { JOB_PATHS } from './contracts';
import type { Job, JobLogPage, JobStatus } from './contracts';
import type { ApiClient, RequestInitLite } from './http';
import type { PageResponse } from '../types/pagination';

/** `active` is queued+running; `terminal` is succeeded/failed/cancelled. */
export type JobStatusFilter = 'active' | 'terminal' | JobStatus;

export interface ListJobsParams {
  status?: JobStatusFilter;
  kind?: string;
  /** 0-based. */
  page?: number;
  /** Serialized as `page_size`. */
  pageSize?: number;
}

/** The seam `useJobEvents` polls through when SSE is unavailable. */
export type JobsFetcher = (
  params: ListJobsParams,
  init?: RequestInitLite,
) => Promise<PageResponse<Job>>;

export const listJobs = (
  c: ApiClient,
  p: ListJobsParams = {},
  init?: RequestInitLite,
): Promise<PageResponse<Job>> =>
  c.get<PageResponse<Job>>(
    JOB_PATHS.list,
    { status: p.status, kind: p.kind, page: p.page, page_size: p.pageSize },
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

/** Stream URL. `since` resumes from a seq, so a reconnect replays missed events. */
export function jobEventsUrl(c: ApiClient, since?: number): string {
  return c.url(JOB_PATHS.events, since && since > 0 ? { since } : undefined);
}
