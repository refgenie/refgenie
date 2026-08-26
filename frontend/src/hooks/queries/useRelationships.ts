import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getRelationshipsExpanded } from '../../services/resources/relationships';

export function useRelationshipsExpanded(digest: string | undefined) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.relationships(digest ?? '', true),
    queryFn: ({ signal }) => getRelationshipsExpanded(client, digest as string, { signal }),
    enabled: !!digest,
  });
}
