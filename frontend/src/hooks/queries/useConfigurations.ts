import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { listConfigurations } from '../../services/resources/configurations';
import type { ListParams } from '../../types/pagination';

export function useConfigurations(p: ListParams = {}, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useResource(
    qk.configurations(p),
    ({ signal }) => listConfigurations(client, p, { signal }),
    { enabled: options?.enabled ?? true },
  );
}
