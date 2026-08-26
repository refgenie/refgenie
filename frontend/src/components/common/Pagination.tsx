import { cn } from '../../utils/cn';
import { currentPage, pageCount } from '../../utils/pagination';
import type { PaginationMeta } from '../../types/pagination';

export interface PaginationProps {
  pagination: PaginationMeta | undefined;
  onOffsetChange: (offset: number) => void;
  /** Numbered pages shown around the current one. */
  window?: number;
}

export function Pagination({ pagination, onOffsetChange, window = 2 }: PaginationProps) {
  if (!pagination || pagination.total === 0) return null;

  const maxPage = pageCount(pagination);
  const page = currentPage(pagination);
  if (maxPage <= 1) return null;

  const first = Math.max(1, page - window);
  const last = Math.min(maxPage, page + window);
  const pages: number[] = [];
  for (let p = first; p <= last; p += 1) pages.push(p);

  const goto = (target: number) => onOffsetChange((target - 1) * pagination.limit);

  return (
    <nav aria-label="Pagination" className="mt-4 flex items-center justify-between gap-4 flex-wrap">
      <p className="text-sm rg-muted">
        {pagination.offset + 1}–{Math.min(pagination.offset + pagination.limit, pagination.total)}{' '}
        of {pagination.total}
      </p>
      <ul className="rg-pagination">
        <li>
          <button
            type="button"
            className="rg-pagination__nav"
            onClick={() => goto(page - 1)}
            disabled={page <= 1}
            aria-label="Previous page"
          >
            ‹
          </button>
        </li>
        {pages.map((p) => (
          <li key={p}>
            <button
              type="button"
              className={cn(
                'rg-pagination__page',
                p === page && 'rg-pagination__page--current',
              )}
              onClick={() => goto(p)}
              aria-current={p === page ? 'page' : undefined}
            >
              {p}
            </button>
          </li>
        ))}
        <li>
          <button
            type="button"
            className="rg-pagination__nav"
            onClick={() => goto(page + 1)}
            disabled={page >= maxPage}
            aria-label="Next page"
          >
            ›
          </button>
        </li>
      </ul>
    </nav>
  );
}
