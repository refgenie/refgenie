import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getRelationshipsExpanded } from '../../services/resources/relationships';

export function useRelationshipsExpanded(digest: string | undefined) {
  const client = useApiClient();
  return useResource(
    qk.relationships(digest ?? '', true),
    ({ signal }) => getRelationshipsExpanded(client, digest as string, { signal }),
    { enabled: !!digest },
  );
}
