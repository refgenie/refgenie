import { describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { Pagination } from './Pagination';
import { buildPageItems, currentPage, pageCount } from '../../utils/pagination';

describe('page math', () => {
  it('rounds up when total is not divisible by limit', () => {
    expect(pageCount({ offset: 0, limit: 50, total: 101 })).toBe(3);
  });

  it('is one page when total fits exactly', () => {
    expect(pageCount({ offset: 0, limit: 50, total: 50 })).toBe(1);
  });

  it('derives the current page from the offset', () => {
    expect(currentPage({ offset: 100, limit: 50, total: 300 })).toBe(3);
  });
});

describe('buildPageItems', () => {
  it('lists every page when they all fit', () => {
    expect(buildPageItems(1, 9)).toEqual([1, 2, 3, 4, 5, 6, 7, 8, 9]);
  });

  it('pins both ends and gaps the middle', () => {
    expect(buildPageItems(8, 15)).toEqual([1, 2, 'gap-left', 7, 8, 9, 'gap-right', 14, 15]);
  });

  it('absorbs the leading gap near the start', () => {
    expect(buildPageItems(1, 15)).toEqual([1, 2, 3, 4, 5, 6, 'gap-right', 14, 15]);
    expect(buildPageItems(5, 15)).toEqual([1, 2, 3, 4, 5, 6, 'gap-right', 14, 15]);
  });

  it('absorbs the trailing gap near the end', () => {
    expect(buildPageItems(15, 15)).toEqual([1, 2, 'gap-left', 10, 11, 12, 13, 14, 15]);
  });

  it('keeps a constant item count across every page of a large set', () => {
    const widths = new Set(
      Array.from({ length: 200 }, (_, i) => buildPageItems(i + 1, 200).length),
    );
    expect([...widths]).toEqual([9]);
  });

  it('never renders a gap that hides only one page', () => {
    for (let maxPage = 10; maxPage <= 40; maxPage += 1) {
      for (let page = 1; page <= maxPage; page += 1) {
        const items = buildPageItems(page, maxPage);
        items.forEach((item, i) => {
          if (typeof item === 'number') return;
          const before = items[i - 1];
          const after = items[i + 1];
          expect(typeof before).toBe('number');
          expect(typeof after).toBe('number');
          // A gap must stand in for at least two pages.
          expect((after as number) - (before as number)).toBeGreaterThan(2);
        });
      }
    }
  });

  it('always shows the current page and never repeats one', () => {
    for (let maxPage = 1; maxPage <= 40; maxPage += 1) {
      for (let page = 1; page <= maxPage; page += 1) {
        const numbers = buildPageItems(page, maxPage).filter(
          (item): item is number => typeof item === 'number',
        );
        expect(numbers).toContain(page);
        expect(numbers).toEqual([...numbers].sort((a, b) => a - b));
        expect(new Set(numbers).size).toBe(numbers.length);
      }
    }
  });

  it('clamps a page hand-edited past the end', () => {
    expect(buildPageItems(99, 15)).toEqual(buildPageItems(15, 15));
    expect(buildPageItems(0, 15)).toEqual(buildPageItems(1, 15));
  });

  it('handles degenerate page counts', () => {
    expect(buildPageItems(1, 1)).toEqual([1]);
    expect(buildPageItems(1, 0)).toEqual([]);
  });

  it('honours custom sibling and boundary counts', () => {
    expect(buildPageItems(10, 20, { siblingCount: 2, boundaryCount: 1 })).toEqual([
      1, 'gap-left', 8, 9, 10, 11, 12, 'gap-right', 20,
    ]);
  });
});

describe('Pagination', () => {
  it('renders nothing when there is a single page', () => {
    const { container } = render(
      <Pagination pagination={{ offset: 0, limit: 50, total: 12 }} onOffsetChange={vi.fn()} />,
    );
    expect(container).toBeEmptyDOMElement();
  });

  it('disables previous on the first page and next on the last', () => {
    const { rerender } = render(
      <Pagination pagination={{ offset: 0, limit: 50, total: 101 }} onOffsetChange={vi.fn()} />,
    );
    expect(screen.getByLabelText('Previous page')).toBeDisabled();
    expect(screen.getByLabelText('Next page')).toBeEnabled();

    rerender(
      <Pagination pagination={{ offset: 100, limit: 50, total: 101 }} onOffsetChange={vi.fn()} />,
    );
    expect(screen.getByLabelText('Previous page')).toBeEnabled();
    expect(screen.getByLabelText('Next page')).toBeDisabled();
  });

  it('reports the offset of the page it navigates to', async () => {
    const onOffsetChange = vi.fn();
    render(
      <Pagination
        pagination={{ offset: 0, limit: 50, total: 300 }}
        onOffsetChange={onOffsetChange}
      />,
    );
    await userEvent.click(screen.getByRole('button', { name: '3' }));
    expect(onOffsetChange).toHaveBeenCalledWith(100);
  });

  it('shows the last page and an ellipsis when there are many', () => {
    const { container } = render(
      <Pagination pagination={{ offset: 0, limit: 10, total: 150 }} onOffsetChange={vi.fn()} />,
    );
    expect(screen.getByRole('button', { name: '15' })).toBeInTheDocument();
    expect(container.querySelectorAll('.rg-pagination__gap')).toHaveLength(1);
    // The ellipsis is decorative, never a target.
    expect(screen.queryByRole('button', { name: '…' })).toBeNull();
  });

  it('jumps to a boundary page at the right offset', async () => {
    const onOffsetChange = vi.fn();
    render(
      <Pagination pagination={{ offset: 0, limit: 10, total: 150 }} onOffsetChange={onOffsetChange} />,
    );
    await userEvent.click(screen.getByRole('button', { name: '15' }));
    expect(onOffsetChange).toHaveBeenCalledWith(140);
  });

  it('does not change width as the user pages through', () => {
    const { container, rerender } = render(
      <Pagination pagination={{ offset: 0, limit: 10, total: 150 }} onOffsetChange={vi.fn()} />,
    );
    const first = container.querySelectorAll('.rg-pagination li').length;
    rerender(
      <Pagination pagination={{ offset: 70, limit: 10, total: 150 }} onOffsetChange={vi.fn()} />,
    );
    expect(container.querySelectorAll('.rg-pagination li').length).toBe(first);
    // 9 page slots plus prev and next.
    expect(first).toBe(11);
  });
});
