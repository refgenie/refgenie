import { walkPages } from '../../utils/pagination';
import { listGenomes } from './genomes';
import type { ApiClient, RequestInitLite } from '../http';
import type { GenomeResponse } from '../../types/api';

/** Every genome present in any of `stores` (server MAX_PAGE_SIZE pages). */
export const fetchCollectionGenomes = (
  c: ApiClient,
  stores: string[],
  init?: RequestInitLite,
): Promise<GenomeResponse[]> =>
  walkPages(
    (page, pageSize) => listGenomes(c, { store: stores, page, pageSize }, init),
    {
      pageSize: 1000,
      maxItems: 20_000,
    },
  );
