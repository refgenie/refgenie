import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { useLocalApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { cancelJob, getJob, getJobLog, listJobs } from '../../services/jobs';
import { useJobStore } from '../../stores/jobStore';
import type { ListJobsParams } from '../../services/jobs';

export function useJobsList(p: ListJobsParams, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useQuery({
    queryKey: qk.jobs(p),
    queryFn: ({ signal }) => listJobs(client, p, { signal }),
    enabled: options?.enabled ?? true,
  });
}

export function useJob(id: string | undefined, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useQuery({
    queryKey: qk.job(id ?? ''),
    queryFn: ({ signal }) => getJob(client, id as string, { signal }),
    enabled: !!id && (options?.enabled ?? true),
  });
}

/** The full log, as opposed to the store's 500-line live ring buffer. */
export function useJobLog(id: string | undefined, offset = 0, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useQuery({
    queryKey: qk.jobLog(id ?? '', offset),
    queryFn: ({ signal }) => getJobLog(client, id as string, offset, { signal }),
    enabled: !!id && (options?.enabled ?? true),
  });
}

export function useCancelJob() {
  const client = useLocalApiClient();
  const queryClient = useQueryClient();
  const upsertJob = useJobStore((state) => state.upsertJob);

  return useMutation({
    mutationFn: (id: string) => cancelJob(client, id),
    onSuccess: (_result, id) => {
      // The authoritative status still arrives on the stream; this only keeps
      // the button from looking inert between the click and the next frame.
      const existing = useJobStore.getState().jobs[id];
      if (existing && existing.status === 'queued') {
        upsertJob({ ...existing, status: 'cancelled' });
      }
      void queryClient.invalidateQueries({ queryKey: ['jobs'] });
    },
  });
}
