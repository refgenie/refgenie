import { useCallback, useMemo, useState } from 'react';
import { useInvalidateResources, useResource } from '../useResource';
import { useLocalApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { cancelJob, getJob, getJobLog, listJobs } from '../../services/jobs';
import { useJobStore } from '../../stores/jobStore';
import type { ListJobsParams } from '../../services/jobs';

export function useJobsList(p: ListJobsParams, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useResource(qk.jobs(p), ({ signal }) => listJobs(client, p, { signal }), {
    enabled: options?.enabled ?? true,
  });
}

export function useJob(id: string | undefined, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useResource(
    qk.job(id ?? ''),
    ({ signal }) => getJob(client, id as string, { signal }),
    { enabled: !!id && (options?.enabled ?? true) },
  );
}

/** The full log, as opposed to the store's 500-line live ring buffer. */
export function useJobLog(id: string | undefined, offset = 0, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useResource(
    qk.jobLog(id ?? '', offset),
    ({ signal }) => getJobLog(client, id as string, offset, { signal }),
    { enabled: !!id && (options?.enabled ?? true) },
  );
}

export interface CancelJobAction {
  mutate: (id: string) => void;
  isPending: boolean;
}

/**
 * The one mutation in the app that is not already a plain async handler.
 * `usePullAction` is the pattern: a `useState` flag and a callback, no
 * mutation framework.
 */
export function useCancelJob(): CancelJobAction {
  const client = useLocalApiClient();
  const invalidate = useInvalidateResources();
  const upsertJob = useJobStore((state) => state.upsertJob);
  const [isPending, setIsPending] = useState(false);

  const mutate = useCallback(
    (id: string) => {
      setIsPending(true);
      cancelJob(client, id)
        .then(() => {
          // The authoritative status still arrives on the stream; this only
          // keeps the button from looking inert between the click and the next
          // frame.
          const existing = useJobStore.getState().jobs[id];
          if (existing && existing.status === 'queued') {
            upsertJob({ ...existing, status: 'cancelled' });
          }
          invalidate(['jobs']);
        })
        // A failed cancel is reported by the job stream, not by this button.
        .catch(() => undefined)
        .finally(() => setIsPending(false));
    },
    [client, invalidate, upsertJob],
  );

  return useMemo(() => ({ mutate, isPending }), [mutate, isPending]);
}
