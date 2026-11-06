import { listQuery } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { AliasPublic } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

export const ALIAS_SEARCH_FIELDS = ['name'] as const;

export interface ListAliasesParams extends ListParams {
  name?: string;
  genome_digest?: string;
}

export const listAliases = (
  c: ApiClient,
  p: ListAliasesParams = {},
  init?: RequestInitLite,
) =>
  c.get<PageResponse<AliasPublic>>(
    '/aliases',
    {
      name: p.name,
      genome_digest: p.genome_digest,
      ...listQuery(p),
    },
    init,
  );
