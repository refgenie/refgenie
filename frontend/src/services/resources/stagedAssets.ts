import { listQuery } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { StagedAssetPublic, StagingMode } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

export interface ListStagedAssetsParams extends ListParams {
  asset_digest?: string;
  mode?: StagingMode;
}

export const listStagedAssets = (
  c: ApiClient,
  p: ListStagedAssetsParams = {},
  init?: RequestInitLite,
) =>
  c.get<PageResponse<StagedAssetPublic>>(
    '/staged_assets',
    {
      asset_digest: p.asset_digest,
      mode: p.mode,
      ...listQuery(p),
    },
    init,
  );
