import { describe, expect, it, vi } from 'vitest';
import { walkPages } from './pagination';
import type { PageResponse } from '../types/pagination';

/** A fake server holding `total` numbered items. */
function server(total: number) {
  return vi.fn(
    async (page: number, pageSize: number): Promise<PageResponse<number>> => {
      const start = page * pageSize;
      const results = Array.from(
        { length: Math.max(0, Math.min(pageSize, total - start)) },
        (_, i) => start + i,
      );
      return { results, pagination: { page, page_size: pageSize, total } };
    },
  );
}

describe('walkPages', () => {
  it('returns nothing for an empty listing after one request', async () => {
    const fetchPage = server(0);
    expect(await walkPages(fetchPage, { pageSize: 10, maxItems: 100 })).toEqual([]);
    expect(fetchPage).toHaveBeenCalledTimes(1);
  });

  it('stops after a single page', async () => {
    const fetchPage = server(3);
    expect(await walkPages(fetchPage, { pageSize: 10, maxItems: 100 })).toEqual([
      0, 1, 2,
    ]);
    expect(fetchPage.mock.calls).toEqual([[0, 10]]);
  });

  it('walks several pages in order', async () => {
    const fetchPage = server(7);
    expect(await walkPages(fetchPage, { pageSize: 3, maxItems: 100 })).toEqual([
      0, 1, 2, 3, 4, 5, 6,
    ]);
    expect(fetchPage.mock.calls).toEqual([
      [0, 3],
      [1, 3],
      [2, 3],
    ]);
  });

  it('stops at the item cap', async () => {
    const fetchPage = server(100);
    expect(await walkPages(fetchPage, { pageSize: 3, maxItems: 5 })).toEqual([
      0, 1, 2, 3, 4, 5,
    ]);
    expect(fetchPage).toHaveBeenCalledTimes(2);
  });
});
