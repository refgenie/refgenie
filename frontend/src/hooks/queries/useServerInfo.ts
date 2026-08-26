import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import {
  getSpeciesSummary,
  getSummary,
  listArchives,
} from '../../services/resources/serverInfo';
import type { ListArchivesParams } from '../../services/resources/serverInfo';

export function useSummary(options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.summary(),
    queryFn: ({ signal }) => getSummary(client, { signal }),
    enabled: options?.enabled ?? true,
    retry: false,
  });
}

export function useSpeciesSummary(options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.speciesSummary(),
    queryFn: ({ signal }) => getSpeciesSummary(client, { signal }),
    enabled: options?.enabled ?? false,
    retry: false,
  });
}

export function useArchives(p: ListArchivesParams, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.archives(p),
    queryFn: ({ signal }) => listArchives(client, p, { signal }),
    enabled: options?.enabled ?? false,
    retry: false,
  });
}
