import type { ApiClient, RequestInitLite } from '../http';
import type { StagedAssetPublic, StagingMode } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

export const STAGED_ASSET_SEARCH_FIELDS = ['asset_digest', 'mode'] as const;

export interface ListStagedAssetsParams extends ListParams {
  asset_digest?: string;
  mode?: StagingMode;
}

export const listStagedAssets = (
  c: ApiClient,
  p: ListStagedAssetsParams = {},
  init?: RequestInitLite,
) =>
  c.get<Paginated<StagedAssetPublic>>(
    '/staged_assets',
    {
      asset_digest: p.asset_digest,
      mode: p.mode,
      q: p.q,
      search_fields: p.searchFields,
      operator: p.operator,
      offset: p.offset,
      limit: p.limit,
    },
    init,
  );

export const getStagedAsset = (c: ApiClient, id: number, init?: RequestInitLite) =>
  c.get<StagedAssetPublic>(`/staged_assets/${id}`, undefined, init);
