/**
 * The app's one job-event stream, mounted by `AppLayout` when `capabilities.jobs`
 * is on. One multiplexed SSE stream serves every job (browsers cap connections
 * per origin). When it drops, the hook re-reads the job list and reopens the
 * stream from the last seq, every few seconds until it holds.
 */

import { useEffect, useRef } from 'react';
import { useLocalApiClient } from './useApiClient';
import { useJobStore } from '../stores/jobStore';
import { jobEventsUrl, listJobs } from '../services/jobs';
import { JOB_EVENT_TYPES } from '../services/contracts';
import type { JobEvent } from '../services/contracts';
import type { JobsFetcher } from '../services/jobs';

/** The subset of `EventSource` this hook uses, so test fakes stay small. */
export interface EventSourceLike {
  addEventListener(type: string, listener: (event: MessageEvent) => void): void;
  close(): void;
  onopen: ((event: Event) => void) | null;
  onerror: ((event: Event) => void) | null;
}

export interface UseJobEventsOptions {
  enabled?: boolean;
  /** Test seam. Defaults to the browser `EventSource`. */
  eventSourceFactory?: (url: string) => EventSourceLike;
  /** Test seam. Defaults to `listJobs` bound to the local client. */
  fetchJobs?: JobsFetcher;
}

/** Matches the "updating every 2s" copy in `JobConsole`. */
const RETRY_MS = 2000;
const SYNC_PAGE_SIZE = 100;

export function useJobEvents({
  enabled = true,
  ...seams
}: UseJobEventsOptions = {}): void {
  const client = useLocalApiClient();
  const seamsRef = useRef(seams);
  seamsRef.current = seams;

  useEffect(() => {
    if (!enabled) return;
    const store = useJobStore.getState;
    const fetchJobs: JobsFetcher =
      seamsRef.current.fetchJobs ?? ((params, init) => listJobs(client, params, init));
    const factory =
      seamsRef.current.eventSourceFactory ?? ((url: string) => new EventSource(url));

    let disposed = false;
    let source: EventSourceLike | null = null;
    let retryTimer: ReturnType<typeof setTimeout> | undefined;

    /** Picks up anything the stream missed while it was down. */
    const sync = async () => {
      try {
        for (const status of ['active', 'terminal'] as const) {
          const page = await fetchJobs({ status, pageSize: SYNC_PAGE_SIZE });
          if (disposed) return;
          store().upsertMany(page.results);
        }
      } catch {
        // The stream, or the next retry, catches up.
      }
    };

    const syncThenOpen = () =>
      void sync().then(() => {
        if (!disposed) open();
      });

    function open() {
      const next = factory(jobEventsUrl(client, store().lastSeq));
      source = next;
      next.onopen = () => store().setConnection('open');
      next.onerror = () => {
        // Reopen by hand rather than let the browser retry, so the URL carries the latest seq.
        next.close();
        store().setConnection('polling');
        retryTimer = setTimeout(syncThenOpen, RETRY_MS);
      };
      for (const type of JOB_EVENT_TYPES) {
        next.addEventListener(type, (event) => {
          let payload: unknown;
          try {
            payload = JSON.parse(String(event.data));
          } catch {
            console.warn('jobs: skipped a malformed SSE frame', event.data);
            return;
          }
          if (payload && typeof payload === 'object') {
            store().applyEvent({ ...payload, type } as JobEvent);
          }
        });
      }
    }

    store().setConnection('connecting');
    syncThenOpen();

    return () => {
      disposed = true;
      clearTimeout(retryTimer);
      source?.close();
    };
  }, [enabled, client]);
}
