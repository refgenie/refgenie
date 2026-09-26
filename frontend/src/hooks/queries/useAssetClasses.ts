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
import type { Paginated } from '../../types/pagination';

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

/** Module-level: passed to `select` by reference so the Maps are built once. */
function toAssetClassIndex(page: Paginated<AssetClassPublic>): AssetClassIndex {
  const byId = new Map<number, AssetClassPublic>();
  const byName = new Map<string, AssetClassPublic>();
  for (const assetClass of page.items) {
    if (assetClass.id !== null && assetClass.id !== undefined) {
      byId.set(assetClass.id, assetClass);
    }
    // AssetClass is unique on (name, version), so a name can repeat across
    // versions. Id lookups are exact; name lookups take the first hit, which is
    // all `input_assets` (which records a name only) can ever resolve to.
    if (!byName.has(assetClass.name)) byName.set(assetClass.name, assetClass);
  }
  return { byId, byName };
}

/**
 * The whole asset class table, indexed both ways. 27 rows on the public server,
 * so one request beats a lookup per reference. `RecipePublic` carries only
 * `output_asset_class_id` and `input_assets` names its inputs by class name, so
 * both directions of the mapping are needed to render and link a recipe.
 */
export function useAssetClassIndex() {
  const client = useApiClient();
  return useResource(
    qk.assetClassIndex(),
    ({ signal }) => listAssetClasses(client, { limit: ASSET_CLASS_INDEX_LIMIT }, { signal }),
    { select: toAssetClassIndex, staleTime: INDEX_STALE_MS },
  );
}
