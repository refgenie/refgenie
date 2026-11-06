import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { NO_RETRY } from '../../services/resourceCache';
import { getAsset, listAssetFiles, listAssets } from '../../services/resources/assets';
import type { ListAssetsParams } from '../../services/resources/assets';

export function useAssets(p: ListAssetsParams, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useResource(qk.assets(p), ({ signal }) => listAssets(client, p, { signal }), {
    enabled: options?.enabled ?? true,
  });
}

export function useAsset(digest: string | undefined) {
  const client = useApiClient();
  return useResource(
    qk.asset(digest ?? ''),
    ({ signal }) => getAsset(client, digest as string, { signal }),
    { enabled: !!digest },
  );
}

export function useAssetFiles(digest: string | undefined, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useResource(
    qk.assetFiles(digest ?? ''),
    ({ signal }) => listAssetFiles(client, digest as string, { signal }),
    { enabled: !!digest && (options?.enabled ?? true), retry: NO_RETRY },
  );
}
