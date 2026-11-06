import { cn } from '../../utils/cn';
import { buildPageItems, pageCount } from '../../utils/pagination';
import type { PageMeta } from '../../types/pagination';

export interface PaginationProps {
  pagination: PageMeta | undefined;
  /** Receives the 0-based page the API takes, not the 1-based label shown. */
  onPageChange: (page: number) => void;
}

export function Pagination({ pagination, onPageChange }: PaginationProps) {
  if (!pagination || pagination.total === 0) return null;

  const maxPage = pageCount(pagination);
  if (maxPage <= 1) return null;

  // Everything rendered is 1-based; only `goto` converts back.
  const page = pagination.page + 1;
  const items = buildPageItems(page, maxPage);
  const goto = (target: number) => onPageChange(target - 1);
  const first = pagination.page * pagination.page_size;

  return (
    <nav aria-label="Pagination" className="mt-4 flex items-center justify-between gap-4 flex-wrap">
      <p className="text-sm rg-muted">
        {first + 1}–{Math.min(first + pagination.page_size, pagination.total)} of{' '}
        {pagination.total}
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
        {items.map((item) =>
          typeof item === 'number' ? (
            <li key={item}>
              <button
                type="button"
                className={cn(
                  'rg-pagination__page',
                  item === page && 'rg-pagination__page--current',
                )}
                onClick={() => goto(item)}
                aria-current={item === page ? 'page' : undefined}
              >
                {item}
              </button>
            </li>
          ) : (
            // Decorative, so screen readers skip it.
            <li key={item} className="rg-pagination__gap" aria-hidden="true">
              …
            </li>
          ),
        )}
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
