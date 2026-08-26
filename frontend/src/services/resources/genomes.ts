import type { ApiClient, RequestInitLite } from '../http';
import type { GenomeDetailResponse, GenomeResponse } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

/** Server-enforced allowlist (`shared.py::list_genomes`); a bad value is a 422. */
export const GENOME_SEARCH_FIELDS = [
  'digest',
  'description',
  'species_name',
  'common_name',
  'assembly_source',
  'assembly_accession',
  'aliases',
] as const;

export interface ListGenomesParams extends ListParams {
  digest?: string;
  alias?: string;
}

export const listGenomes = (
  c: ApiClient,
  p: ListGenomesParams = {},
  init?: RequestInitLite,
) =>
  c.get<Paginated<GenomeResponse>>(
    '/genomes',
    {
      digest: p.digest,
      alias: p.alias,
      q: p.q,
      search_fields: p.searchFields,
      operator: p.operator,
      offset: p.offset,
      limit: p.limit,
    },
    init,
  );

export const getGenome = (c: ApiClient, digest: string, init?: RequestInitLite) =>
  c.get<GenomeDetailResponse>(
    `/genomes/${encodeURIComponent(digest)}`,
    undefined,
    init,
  );
