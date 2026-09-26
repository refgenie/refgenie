import type { ApiClient, RequestInitLite } from '../http';
import type { GenomeDetailResponse, GenomeResponse } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

/** Server-enforced allowlist (`catalog.py::list_genomes`); a bad value is a 422. */
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

/** The server's MAX_PAGE_SIZE (`refgenie/const.py`). */
const INDEX_PAGE_SIZE = 1000;

/**
 * Hard stop for the bulk index. A backend past this size needs a server-side
 * alias join (see the NOT-APPROVED section of the class-first browsing plan),
 * not a longer loop.
 */
const INDEX_MAX_GENOMES = 20_000;

/**
 * Every genome on this backend, in one pass.
 *
 * `/asset_groups` and `/assets` carry only `genome_digest`, and there is no bulk
 * "resolve these digests" endpoint — `/genomes?digest=` and
 * `/aliases?genome_digest=` each take exactly one. Fetching the whole list once
 * (701 rows / ~55 KB gzipped on the public server) is a single request; the
 * alternative is one request per row on every page flip.
 */
export async function fetchAllGenomes(
  c: ApiClient,
  init?: RequestInitLite,
): Promise<GenomeResponse[]> {
  const all: GenomeResponse[] = [];
  let offset = 0;
  for (;;) {
    const page = await listGenomes(c, { offset, limit: INDEX_PAGE_SIZE }, init);
    all.push(...page.items);
    if (page.items.length === 0) break;
    offset += INDEX_PAGE_SIZE;
    if (offset >= page.pagination.total) break;
    if (all.length >= INDEX_MAX_GENOMES) break;
  }
  return all;
}
