import { listQuery } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { AssetClassPublic, AssetClassResponse } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

export const ASSET_CLASS_SEARCH_FIELDS = ['name', 'version', 'description'] as const;

export interface ListAssetClassesParams extends ListParams {
  name?: string;
  version?: string;
}

export const listAssetClasses = (
  c: ApiClient,
  p: ListAssetClassesParams = {},
  init?: RequestInitLite,
) =>
  c.get<PageResponse<AssetClassPublic>>(
    '/asset_classes',
    {
      name: p.name,
      version: p.version,
      ...listQuery(p),
    },
    init,
  );

export const getAssetClass = (c: ApiClient, id: number, init?: RequestInitLite) =>
  c.get<AssetClassResponse>(`/asset_classes/${id}`, undefined, init);

/** The whole asset class table in one request; the server caps at 1000. */
export const ASSET_CLASS_INDEX_LIMIT = 1000;
