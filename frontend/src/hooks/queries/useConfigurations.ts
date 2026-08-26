import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { listConfigurations } from '../../services/resources/configurations';
import type { ListParams } from '../../types/pagination';

export function useConfigurations(p: ListParams = {}, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.configurations(p),
    queryFn: ({ signal }) => listConfigurations(client, p, { signal }),
    enabled: options?.enabled ?? true,
  });
}
