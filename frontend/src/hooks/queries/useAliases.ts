import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getAlias, listAliases } from '../../services/resources/aliases';
import type { ListAliasesParams } from '../../services/resources/aliases';

export function useAliases(p: ListAliasesParams) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.aliases(p),
    queryFn: ({ signal }) => listAliases(client, p, { signal }),
  });
}

export function useAlias(name: string | undefined) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.alias(name ?? ''),
    queryFn: ({ signal }) => getAlias(client, name as string, { signal }),
    enabled: !!name,
  });
}
