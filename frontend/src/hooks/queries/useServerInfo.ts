import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { NO_RETRY } from '../../services/resourceCache';
import { listArchives } from '../../services/resources/serverInfo';
import type { ListArchivesParams } from '../../services/resources/serverInfo';

export function useArchives(p: ListArchivesParams, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useResource(qk.archives(p), ({ signal }) => listArchives(client, p, { signal }), {
    enabled: options?.enabled ?? false,
    retry: NO_RETRY,
  });
}
