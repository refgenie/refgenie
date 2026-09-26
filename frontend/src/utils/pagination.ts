import type { PaginationMeta } from '../types/pagination';

export function pageCount(meta: PaginationMeta): number {
  if (meta.limit <= 0) return 1;
  return Math.max(1, Math.ceil(meta.total / meta.limit));
}

export function currentPage(meta: PaginationMeta): number {
  if (meta.limit <= 0) return 1;
  return Math.floor(meta.offset / meta.limit) + 1;
}

/**
 * One slot in the rendered page list: either a page number, or a gap standing
 * in for two or more hidden pages. There is at most one gap on each side, so
 * the two string literals double as stable React keys.
 */
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
 * Build the page list for a paginator: pinned pages at both ends, a run that
 * slides with the current page, and gaps between them — `1 2 … 7 8 9 … 14 15`.
 *
 * Two properties the tests pin down, because both are what make the control
 * usable rather than merely correct:
 *
 * 1. **Constant width.** Whenever `maxPage` is large enough to need a gap, the
 *    result is always `boundaryCount * 2 + siblingCount * 2 + 3` items. Buttons
 *    do not move under the cursor as the user pages through.
 * 2. **No gap ever hides a single page.** A `…` standing in for one page is
 *    worse than the page itself, so in that case the run absorbs it.
 */
export function buildPageItems(
  page: number,
  maxPage: number,
  { siblingCount = 1, boundaryCount = 2 }: PageItemOptions = {},
): PageItem[] {
  // The sliding run: the current page plus its siblings.
  const runLength = siblingCount * 2 + 1;
  // Two boundaries, two gaps, and the run.
  const totalSlots = boundaryCount * 2 + runLength + 2;

  // Few enough pages that every one fits. No gaps, nothing hidden.
  if (maxPage <= totalSlots) return range(1, maxPage);

  // An offset hand-edited past the end must not produce a broken control.
  const active = Math.min(Math.max(page, 1), maxPage);

  const firstInterior = boundaryCount + 1;
  const lastInterior = maxPage - boundaryCount;

  let start = active - siblingCount;
  let end = active + siblingCount;

  // When the run reaches a boundary, the gap on that side would hide zero or
  // one page. Drop it and let the run grow by exactly one, which keeps the
  // item count at `totalSlots`. Only one side can trigger: both would require
  // `maxPage <= totalSlots`, and that returned above.
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
