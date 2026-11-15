import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { listAliases } from '../../services/resources/aliases';
import type { ListAliasesParams } from '../../services/resources/aliases';

export function useAliases(p: ListAliasesParams, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useResource(
    qk.aliases(p),
    ({ signal }) => listAliases(client, p, { signal }),
    { enabled: options?.enabled ?? true },
  );
}
