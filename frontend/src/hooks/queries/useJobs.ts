import { useCallback, useMemo, useState } from 'react';
import { useInvalidate, useResource } from '../useResource';
import { useLocalApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { cancelJob, getJobLog, listJobs } from '../../services/jobs';
import { useJobStore } from '../../stores/jobStore';
import type { ListJobsParams } from '../../services/jobs';

export function useJobsList(p: ListJobsParams, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useResource(qk.jobs(p), ({ signal }) => listJobs(client, p, { signal }), {
    enabled: options?.enabled ?? true,
  });
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

/** Cancels a job, with a pending flag for the button. */
export function useCancelJob(): CancelJobAction {
  const client = useLocalApiClient();
  const invalidate = useInvalidate();
  const upsertJob = useJobStore((state) => state.upsertJob);
  const [isPending, setIsPending] = useState(false);

  const mutate = useCallback(
    (id: string) => {
      setIsPending(true);
      cancelJob(client, id)
        .then(() => {
          // The stream sends the real status; this just updates the button sooner.
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
