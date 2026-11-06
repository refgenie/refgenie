import { listQuery } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { AssetFilesResponse, AssetResponse } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

/**
 * At most one of these filters: the server joins AssetGroup once per filter, so
 * any two give a 500. Flatten this once the server is fixed.
 */
type AssetGroupFilter =
  | { genome_digest?: string; asset_group_name?: never; asset_group_id?: never }
  | { genome_digest?: never; asset_group_name?: string; asset_group_id?: never }
  | { genome_digest?: never; asset_group_name?: never; asset_group_id?: number };

/**
 * `q` here searches only the asset tag (`default`), not the group name users
 * type (`fasta`), so the UI searches asset classes and groups instead.
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
  c.get<PageResponse<AssetResponse>>(
    '/assets',
    {
      name: p.name,
      asset_group_name: p.asset_group_name,
      genome_digest: p.genome_digest,
      recipe_name: p.recipe_name,
      asset_group_id: p.asset_group_id,
      ...listQuery(p),
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
