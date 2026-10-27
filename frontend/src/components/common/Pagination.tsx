import { cn } from '../../utils/cn';
import { buildPageItems, currentPage, pageCount } from '../../utils/pagination';
import type { PaginationMeta } from '../../types/pagination';

export interface PaginationProps {
  pagination: PaginationMeta | undefined;
  onOffsetChange: (offset: number) => void;
}

export function Pagination({ pagination, onOffsetChange }: PaginationProps) {
  if (!pagination || pagination.total === 0) return null;

  const maxPage = pageCount(pagination);
  const page = currentPage(pagination);
  if (maxPage <= 1) return null;

  const items = buildPageItems(page, maxPage);
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
            // Decorative: the pages it stands for are not reachable from here,
            // so announcing "horizontal ellipsis" is noise.
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
