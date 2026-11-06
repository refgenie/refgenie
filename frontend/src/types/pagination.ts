/** Pagination and search wire types (`refgenie/server/pagination.py`, `search.py`). */

/** `page` is 0-based; the server returns 400 for a page past the end. */
export interface PageMeta {
  page: number;
  page_size: number;
  total: number;
}

export interface PageResponse<T> {
  results: T[];
  pagination: PageMeta;
}

/** `eq` is case-sensitive; the others are not. The server default is `contains`. */
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
  /** 0-based. Default 0. */
  page?: number;
  /** Serialized as `page_size`. Server default 100, max 1000. */
  pageSize?: number;
}

/** UI page size. The server cap is MAX_PAGE_SIZE = 1000. */
export const DEFAULT_PAGE_SIZE = 50;
