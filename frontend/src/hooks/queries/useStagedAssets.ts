import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { listStagedAssets } from '../../services/resources/stagedAssets';
import type { ListStagedAssetsParams } from '../../services/resources/stagedAssets';

export function useStagedAssets(
  p: ListStagedAssetsParams,
  options?: { enabled?: boolean },
) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.stagedAssets(p),
    queryFn: ({ signal }) => listStagedAssets(client, p, { signal }),
    enabled: options?.enabled ?? true,
  });
}
