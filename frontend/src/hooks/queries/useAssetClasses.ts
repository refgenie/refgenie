import { useMemo } from 'react';
import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import {
  ASSET_CLASS_INDEX_LIMIT,
  getAssetClass,
  listAssetClasses,
} from '../../services/resources/assetClasses';
import type { ListAssetClassesParams } from '../../services/resources/assetClasses';
import type { AssetClassPublic } from '../../types/api';
import type { PageResponse } from '../../types/pagination';

export function useAssetClasses(p: ListAssetClassesParams) {
  const client = useApiClient();
  return useResource(qk.assetClasses(p), ({ signal }) =>
    listAssetClasses(client, p, { signal }),
  );
}

export function useAssetClass(id: number | undefined) {
  const client = useApiClient();
  return useResource(
    qk.assetClass(id ?? -1),
    ({ signal }) => getAssetClass(client, id as number, { signal }),
    { enabled: id !== undefined && Number.isFinite(id) },
  );
}

const INDEX_STALE_MS = 5 * 60 * 1000;

export interface AssetClassIndex {
  byId: Map<number, AssetClassPublic>;
  byName: Map<string, AssetClassPublic>;
}

function toAssetClassIndex(page: PageResponse<AssetClassPublic>): AssetClassIndex {
  const byId = new Map<number, AssetClassPublic>();
  const byName = new Map<string, AssetClassPublic>();
  for (const assetClass of page.results) {
    if (assetClass.id !== null && assetClass.id !== undefined) {
      byId.set(assetClass.id, assetClass);
    }
    // A name can repeat across versions; name lookups take the first hit.
    if (!byName.has(assetClass.name)) byName.set(assetClass.name, assetClass);
  }
  return { byId, byName };
}

/**
 * The whole (small) asset class table, indexed by id and by name, since recipes
 * name outputs by id and inputs by name.
 */
export function useAssetClassIndex() {
  const client = useApiClient();
  const query = useResource(
    qk.assetClassIndex(),
    ({ signal }) => listAssetClasses(client, { pageSize: ASSET_CLASS_INDEX_LIMIT }, { signal }),
    { staleTime: INDEX_STALE_MS },
  );
  const data = useMemo(() => query.data && toAssetClassIndex(query.data), [query.data]);
  return { ...query, data };
}
