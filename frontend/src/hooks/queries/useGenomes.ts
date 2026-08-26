import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getGenome, listGenomes } from '../../services/resources/genomes';
import type { ListGenomesParams } from '../../services/resources/genomes';

export function useGenomes(p: ListGenomesParams) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.genomes(p),
    queryFn: ({ signal }) => listGenomes(client, p, { signal }),
  });
}

export function useGenome(digest: string | undefined) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.genome(digest ?? ''),
    queryFn: ({ signal }) => getGenome(client, digest as string, { signal }),
    enabled: !!digest,
  });
}
