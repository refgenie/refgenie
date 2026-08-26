/**
 * Pagination and search wire types.
 *
 * Mirrors `refgenie/utils/pagination.py` (PaginationMeta / PaginatedResponse)
 * and `refgenie/utils/search.py` (SearchOperator / SearchParams).
 */

export interface PaginationMeta {
  offset: number;
  limit: number;
  total: number;
}

export interface Paginated<T> {
  items: T[];
  pagination: PaginationMeta;
}

/**
 * `eq` is a case-sensitive `==`; the other three are case-insensitive SQL
 * `ilike`. The server default is `contains`.
 */
export type SearchOperator = 'eq' | 'contains' | 'starts_with' | 'ends_with';

export const SEARCH_OPERATORS: readonly SearchOperator[] = [
  'contains',
  'eq',
  'starts_with',
  'ends_with',
] as const;

export interface ListParams {
  q?: string;
  /** Serialized as a comma-joined `search_fields` query parameter. */
  searchFields?: string[];
  operator?: SearchOperator;
  /** Default 0. */
  offset?: number;
  /** Server default 100, server cap 1000 (refgenie/const.py). */
  limit?: number;
}

/** UI page size. The server cap is MAX_PAGE_SIZE = 1000. */
export const DEFAULT_PAGE_SIZE = 50;
