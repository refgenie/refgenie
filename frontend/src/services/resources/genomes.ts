import { listQuery, walkPages } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { GenomeDetailResponse, GenomeResponse } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

/** Fields the server allows for search (`catalog.py::list_genomes`). */
export const GENOME_SEARCH_FIELDS = [
  'digest',
  'description',
  'species',
  'synonyms',
  'assembly_source',
  'assembly_accession',
  'aliases',
] as const;

export interface ListGenomesParams extends ListParams {
  digest?: string;
  alias?: string;
  has_assets?: boolean;
  taxon_id?: number;
}

export const listGenomes = (
  c: ApiClient,
  p: ListGenomesParams = {},
  init?: RequestInitLite,
) =>
  c.get<PageResponse<GenomeResponse>>(
    '/genomes',
    {
      digest: p.digest,
      alias: p.alias,
      has_assets: p.has_assets,
      taxon_id: p.taxon_id,
      ...listQuery(p),
    },
    init,
  );

export const getGenome = (c: ApiClient, digest: string, init?: RequestInitLite) =>
  c.get<GenomeDetailResponse>(
    `/genomes/${encodeURIComponent(digest)}`,
    undefined,
    init,
  );

/** The server's MAX_PAGE_SIZE (`refgenie/const.py`). */
const INDEX_PAGE_SIZE = 1000;

/** Cap on the bulk index; a larger backend needs a server-side join instead. */
const INDEX_MAX_GENOMES = 20_000;

/**
 * Fetch every genome. Asset lists carry only `genome_digest` and there is no
 * bulk lookup endpoint, so one full fetch beats a request per row.
 */
export const fetchAllGenomes = (
  c: ApiClient,
  init?: RequestInitLite,
): Promise<GenomeResponse[]> =>
  walkPages((page, pageSize) => listGenomes(c, { page, pageSize }, init), {
    pageSize: INDEX_PAGE_SIZE,
    maxItems: INDEX_MAX_GENOMES,
  });
