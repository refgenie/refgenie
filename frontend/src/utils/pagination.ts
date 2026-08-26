import type { PaginationMeta } from '../types/pagination';

export function pageCount(meta: PaginationMeta): number {
  if (meta.limit <= 0) return 1;
  return Math.max(1, Math.ceil(meta.total / meta.limit));
}

export function currentPage(meta: PaginationMeta): number {
  if (meta.limit <= 0) return 1;
  return Math.floor(meta.offset / meta.limit) + 1;
}
