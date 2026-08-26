import type { ApiClient, RequestInitLite } from '../http';
import type { AssetFilesResponse, AssetResponse } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

export const ASSET_SEARCH_FIELDS = ['name', 'digest', 'path'] as const;

export interface ListAssetsParams extends ListParams {
  name?: string;
  asset_group_name?: string;
  genome_digest?: string;
  recipe_name?: string;
  asset_group_id?: number;
}

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
