import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getAssetGroup, listAssetGroups } from '../../services/resources/assetGroups';
import type { ListAssetGroupsParams } from '../../services/resources/assetGroups';

export function useAssetGroups(p: ListAssetGroupsParams, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.assetGroups(p),
    queryFn: ({ signal }) => listAssetGroups(client, p, { signal }),
    enabled: options?.enabled ?? true,
  });
}

export function useAssetGroup(id: number | undefined) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.assetGroup(id ?? -1),
    queryFn: ({ signal }) => getAssetGroup(client, id as number, { signal }),
    enabled: id !== undefined && Number.isFinite(id),
  });
}
