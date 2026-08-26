import type { ApiClient, RequestInitLite } from '../http';
import type { AssetClassPublic } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

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
  c.get<Paginated<AssetClassPublic>>(
    '/asset_classes',
    {
      name: p.name,
      version: p.version,
      q: p.q,
      search_fields: p.searchFields,
      operator: p.operator,
      offset: p.offset,
      limit: p.limit,
    },
    init,
  );

export const getAssetClass = (c: ApiClient, id: number, init?: RequestInitLite) =>
  c.get<AssetClassPublic>(`/asset_classes/${id}`, undefined, init);
