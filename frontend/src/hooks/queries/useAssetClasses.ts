import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getAssetClass, listAssetClasses } from '../../services/resources/assetClasses';
import type { ListAssetClassesParams } from '../../services/resources/assetClasses';

export function useAssetClasses(p: ListAssetClassesParams) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.assetClasses(p),
    queryFn: ({ signal }) => listAssetClasses(client, p, { signal }),
  });
}

export function useAssetClass(id: number | undefined) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.assetClass(id ?? -1),
    queryFn: ({ signal }) => getAssetClass(client, id as number, { signal }),
    enabled: id !== undefined && Number.isFinite(id),
  });
}
