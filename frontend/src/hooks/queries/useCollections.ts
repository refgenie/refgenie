import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { useUiConfig } from '../useUiConfig';
import { useGenomes } from './useGenomes';
import { qk } from '../../services/queryKeys';
import { fetchCollectionGenomes } from '../../services/resources/collections';
import type { CollectionConfig } from '../../types/ui';

const COLLECTION_STALE_MS = 5 * 60 * 1000;

/** Collections the server declares in `/service-info`. */
export function useCollections(): CollectionConfig[] {
  return useUiConfig().collections;
}

export function useCollection(id: string | undefined): CollectionConfig | undefined {
  const collections = useCollections();
  return id ? collections.find((c) => c.id === id) : undefined;
}

/**
 * Every genome in the collection, as one cache entry. Keyed in the `genomes`
 * family so pull/build/delete invalidation refreshes it too.
 */
export function useCollectionGenomes(collection: CollectionConfig | undefined) {
  const client = useApiClient();
  const stores = collection?.stores ?? [];
  return useResource(
    qk.genomes({ collection: collection?.id, stores, all: true }),
    ({ signal }) => fetchCollectionGenomes(client, stores, { signal }),
    { enabled: !!collection, staleTime: COLLECTION_STALE_MS },
  );
}

/** Genome count for an index card: one row requested, `pagination.total` read. */
export function useCollectionCount(collection: CollectionConfig) {
  const query = useGenomes({ store: collection.stores, pageSize: 1 });
  return { ...query, data: query.data?.pagination.total };
}
