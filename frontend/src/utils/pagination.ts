import type { ListParams, PageMeta, PageResponse } from '../types/pagination';

export function pageCount(meta: PageMeta): number {
  if (meta.page_size <= 0) return 1;
  return Math.max(1, Math.ceil(meta.total / meta.page_size));
}

/** True when no page follows this one; the stop condition for a fetch-all loop. */
export function isLastPage(meta: PageMeta): boolean {
  return (meta.page + 1) * meta.page_size >= meta.total;
}

export interface WalkPagesOptions {
  pageSize: number;
  /** Stop once this many items are collected; the last page is kept whole. */
  maxItems: number;
}

/** Fetch pages from 0 until the last page or the item cap; return every item. */
export async function walkPages<T>(
  fetchPage: (page: number, pageSize: number) => Promise<PageResponse<T>>,
  { pageSize, maxItems }: WalkPagesOptions,
): Promise<T[]> {
  const all: T[] = [];
  for (let page = 0; ; page += 1) {
    const data = await fetchPage(page, pageSize);
    all.push(...data.results);
    if (isLastPage(data.pagination) || all.length >= maxItems) return all;
  }
}

/** A page number, or a gap hiding two or more pages (also used as React keys). */
export type PageItem = number | 'gap-left' | 'gap-right';

export interface PageItemOptions {
  /** Pages shown either side of the current one. */
  siblingCount?: number;
  /** Pages pinned at each end. */
  boundaryCount?: number;
}

function range(from: number, to: number): number[] {
  const out: number[] = [];
  for (let n = from; n <= to; n += 1) out.push(n);
  return out;
}

/**
 * Page list for a paginator, e.g. `1 2 … 7 8 9 … 14 15`. The item count stays
 * constant once gaps appear, so buttons do not shift, and a gap never hides just
 * one page.
 */
export function buildPageItems(
  page: number,
  maxPage: number,
  { siblingCount = 1, boundaryCount = 2 }: PageItemOptions = {},
): PageItem[] {
  const runLength = siblingCount * 2 + 1;
  // Two boundaries, two gaps, and the run.
  const totalSlots = boundaryCount * 2 + runLength + 2;

  if (maxPage <= totalSlots) return range(1, maxPage);

  // Clamp a page typed past the end.
  const active = Math.min(Math.max(page, 1), maxPage);

  const firstInterior = boundaryCount + 1;
  const lastInterior = maxPage - boundaryCount;

  let start = active - siblingCount;
  let end = active + siblingCount;

  // Near a boundary the gap would hide zero or one page: drop it and grow the run
  // by one to keep the count at `totalSlots`.
  if (start <= firstInterior + 1) {
    start = firstInterior;
    end = firstInterior + runLength;
  } else if (end >= lastInterior - 1) {
    end = lastInterior;
    start = lastInterior - runLength;
  }

  return [
    ...range(1, boundaryCount),
    ...(start > firstInterior ? (['gap-left'] as PageItem[]) : []),
    ...range(start, end),
    ...(end < lastInterior ? (['gap-right'] as PageItem[]) : []),
    ...range(maxPage - boundaryCount + 1, maxPage),
  ];
}

/** Wire names for the shared list and search parameters. */
export function listQuery(p: ListParams) {
  return {
    q: p.q,
    search_fields: p.searchFields,
    operator: p.operator,
    page: p.page,
    page_size: p.pageSize,
  };
}
