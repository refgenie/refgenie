import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getAlias, listAliases } from '../../services/resources/aliases';
import type { ListAliasesParams } from '../../services/resources/aliases';

export function useAliases(p: ListAliasesParams) {
  const client = useApiClient();
  return useResource(qk.aliases(p), ({ signal }) => listAliases(client, p, { signal }));
}

export function useAlias(name: string | undefined) {
  const client = useApiClient();
  return useResource(
    qk.alias(name ?? ''),
    ({ signal }) => getAlias(client, name as string, { signal }),
    { enabled: !!name },
  );
}
