import { useResource } from '../useResource';
import { useLocalApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { NO_RETRY } from '../../services/resourceCache';
import { getLocalFilePreview, listLocalAssetFiles } from '../../services/resources/localFiles';

export function useLocalAssetFiles(digest: string | undefined, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useResource(
    qk.localAssetFiles(digest ?? ''),
    ({ signal }) => listLocalAssetFiles(client, digest as string, { signal }),
    { enabled: !!digest && (options?.enabled ?? true), retry: NO_RETRY },
  );
}

export function useLocalFilePreview(
  digest: string,
  path: string | null,
  options?: { enabled?: boolean },
) {
  const client = useLocalApiClient();
  return useResource(
    qk.localFilePreview(digest, path ?? ''),
    ({ signal }) => getLocalFilePreview(client, digest, path as string, { signal }),
    { enabled: !!digest && !!path && (options?.enabled ?? true), retry: NO_RETRY },
  );
}
