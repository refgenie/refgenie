import { describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { Pagination } from './Pagination';
import { currentPage, pageCount } from '../../utils/pagination';

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
});
