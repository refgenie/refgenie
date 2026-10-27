import type { ApiClient, RequestInitLite } from '../http';
import type { AliasPublic, AliasResponse } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

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
  c.get<Paginated<AliasPublic>>(
    '/aliases',
    {
      name: p.name,
      genome_digest: p.genome_digest,
      q: p.q,
      search_fields: p.searchFields,
      operator: p.operator,
      offset: p.offset,
      limit: p.limit,
    },
    init,
  );

export const getAlias = (c: ApiClient, name: string, init?: RequestInitLite) =>
  c.get<AliasResponse>(`/aliases/${encodeURIComponent(name)}`, undefined, init);
