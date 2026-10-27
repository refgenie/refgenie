import type { ApiClient, RequestInitLite } from '../http';
import type { AssetFilesResponse, AssetResponse } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

/**
 * `catalog.py::list_assets` calls `query.join(AssetGroup)` once per filter, so
 * ANY TWO of `genome_digest`, `asset_group_name`, `asset_group_id` produce a
 * duplicate join and a 500. Verified live against api.refgenie.org: each alone
 * is a 200, every pair is a 500. `name` and `recipe_name` join other tables and
 * combine freely.
 *
 * The union makes the failing combination unrepresentable rather than merely
 * undocumented. Collapse it back into a flat interface once the server joins
 * AssetGroup at most once.
 */
type AssetGroupFilter =
  | { genome_digest?: string; asset_group_name?: never; asset_group_id?: never }
  | { genome_digest?: never; asset_group_name?: string; asset_group_id?: never }
  | { genome_digest?: never; asset_group_name?: never; asset_group_id?: number };

/**
 * `/v4/assets` searches `name`, `digest` and `path` only, and its `name` is the
 * TAG (`default`, `2.5.5`) — the word a human types (`fasta`, `bowtie2_index`)
 * lives in `asset_group_name`, which the endpoint cannot search. Free-text
 * asset search therefore goes through `/v4/asset_classes` and
 * `/v4/asset_groups`; no UI wires `q` to this endpoint.
 */
export type ListAssetsParams = ListParams & {
  name?: string;
  recipe_name?: string;
} & AssetGroupFilter;

export const listAssets = (
  c: ApiClient,
  p: ListAssetsParams = {},
  init?: RequestInitLite,
) =>
  c.get<Paginated<AssetResponse>>(
    '/assets',
    {
      name: p.name,
      asset_group_name: p.asset_group_name,
      genome_digest: p.genome_digest,
      recipe_name: p.recipe_name,
      asset_group_id: p.asset_group_id,
      q: p.q,
      search_fields: p.searchFields,
      operator: p.operator,
      offset: p.offset,
      limit: p.limit,
    },
    init,
  );

export const getAsset = (c: ApiClient, digest: string, init?: RequestInitLite) =>
  c.get<AssetResponse>(`/assets/${encodeURIComponent(digest)}`, undefined, init);

export const listAssetFiles = (c: ApiClient, digest: string, init?: RequestInitLite) =>
  c.get<AssetFilesResponse>(
    `/assets/${encodeURIComponent(digest)}/files`,
    undefined,
    init,
  );
