import { listQuery } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { AssetGroupPublic } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

export interface ListAssetGroupsParams extends ListParams {
  genome_digest?: string;
  asset_class?: string;
  asset_group_name?: string;
  asset_group_id?: number;
}

export const listAssetGroups = (
  c: ApiClient,
  p: ListAssetGroupsParams = {},
  init?: RequestInitLite,
) =>
  c.get<PageResponse<AssetGroupPublic>>(
    '/asset_groups',
    {
      genome_digest: p.genome_digest,
      asset_class: p.asset_class,
      asset_group_name: p.asset_group_name,
      asset_group_id: p.asset_group_id,
      ...listQuery(p),
    },
    init,
  );

export const getAssetGroup = (c: ApiClient, id: number, init?: RequestInitLite) =>
  c.get<AssetGroupPublic>(`/asset_groups/${id}`, undefined, init);
