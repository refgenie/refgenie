import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getAsset, listAssetFiles, listAssets } from '../../services/resources/assets';
import type { ListAssetsParams } from '../../services/resources/assets';

export function useAssets(p: ListAssetsParams, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.assets(p),
    queryFn: ({ signal }) => listAssets(client, p, { signal }),
    enabled: options?.enabled ?? true,
  });
}

export function useAsset(digest: string | undefined) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.asset(digest ?? ''),
    queryFn: ({ signal }) => getAsset(client, digest as string, { signal }),
    enabled: !!digest,
  });
}

export function useAssetFiles(digest: string | undefined, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.assetFiles(digest ?? ''),
    queryFn: ({ signal }) => listAssetFiles(client, digest as string, { signal }),
    enabled: !!digest && (options?.enabled ?? true),
    retry: false,
  });
}
