import { describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { Pagination } from './Pagination';
import { buildPageItems, isLastPage, pageCount } from '../../utils/pagination';
import type { PageItemOptions } from '../../utils/pagination';

describe('page math', () => {
  it.each([
    [101, 3, [false, false, true]],
    [100, 2, [false, true]],
    [50, 1, [true]],
    [0, 1, [true]],
  ])('total %i at page size 50: %i page(s)', (total, count, last) => {
    expect(pageCount({ page: 0, page_size: 50, total })).toBe(count);
    expect(last.map((_, page) => isLastPage({ page, page_size: 50, total }))).toEqual(
      last,
    );
  });
});

describe('buildPageItems', () => {
  it.each<[number, number, PageItemOptions, ReturnType<typeof buildPageItems>]>([
    [1, 9, {}, [1, 2, 3, 4, 5, 6, 7, 8, 9]],
    [8, 15, {}, [1, 2, 'gap-left', 7, 8, 9, 'gap-right', 14, 15]],
    [1, 15, {}, [1, 2, 3, 4, 5, 6, 'gap-right', 14, 15]],
    [5, 15, {}, [1, 2, 3, 4, 5, 6, 'gap-right', 14, 15]],
    [15, 15, {}, [1, 2, 'gap-left', 10, 11, 12, 13, 14, 15]],
    [99, 15, {}, [1, 2, 'gap-left', 10, 11, 12, 13, 14, 15]],
    [0, 15, {}, [1, 2, 3, 4, 5, 6, 'gap-right', 14, 15]],
    [1, 1, {}, [1]],
    [1, 0, {}, []],
    [
      10,
      20,
      { siblingCount: 2, boundaryCount: 1 },
      [1, 'gap-left', 8, 9, 10, 11, 12, 'gap-right', 20],
    ],
  ])('page %i of %i %o', (page, maxPage, options, expected) => {
    expect(buildPageItems(page, maxPage, options)).toEqual(expected);
  });

  it('keeps a constant width, shows the current page, and never gaps a single page', () => {
    for (let maxPage = 10; maxPage <= 40; maxPage += 1) {
      for (let page = 1; page <= maxPage; page += 1) {
        const items = buildPageItems(page, maxPage);
        const numbers = items.filter(
          (item): item is number => typeof item === 'number',
        );
        expect(items).toHaveLength(9);
        expect(numbers).toContain(page);
        expect(numbers).toEqual([...new Set(numbers)].sort((a, b) => a - b));
        items.forEach((item, i) => {
          if (typeof item !== 'number') {
            expect((items[i + 1] as number) - (items[i - 1] as number)).toBeGreaterThan(
              2,
            );
          }
        });
      }
    }
  });
});

describe('Pagination', () => {
  it('renders nothing when there is a single page', () => {
    const { container } = render(
      <Pagination
        pagination={{ page: 0, page_size: 50, total: 12 }}
        onPageChange={vi.fn()}
      />,
    );
    expect(container).toBeEmptyDOMElement();
  });

  it('shows 1-based pages but reports the 0-based page the API takes', async () => {
    const onPageChange = vi.fn();
    const { container, rerender } = render(
      <Pagination
        pagination={{ page: 0, page_size: 10, total: 150 }}
        onPageChange={onPageChange}
      />,
    );
    expect(screen.getByRole('button', { name: '1' })).toHaveAttribute(
      'aria-current',
      'page',
    );
    expect(screen.getByText(/1–10 of/)).toBeInTheDocument();
    expect(screen.getByLabelText('Previous page')).toBeDisabled();
    expect(container.querySelectorAll('.rg-pagination__gap')).toHaveLength(1);

    await userEvent.click(screen.getByRole('button', { name: '15' }));
    expect(onPageChange).toHaveBeenLastCalledWith(14);

    rerender(
      <Pagination
        pagination={{ page: 14, page_size: 10, total: 150 }}
        onPageChange={onPageChange}
      />,
    );
    expect(screen.getByText(/141–150 of/)).toBeInTheDocument();
    expect(screen.getByLabelText('Next page')).toBeDisabled();
    await userEvent.click(screen.getByLabelText('Previous page'));
    expect(onPageChange).toHaveBeenLastCalledWith(13);
  });
});
