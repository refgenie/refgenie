/**
 * Returns to the first page when the server answers 400 for a page past the end
 * (a stale link, or deleting the last row). A 400 on page 0 is a real error.
 */

import { useEffect } from 'react';
import { ApiError } from '../services/http';

export function useResetPastEndPage(
  error: unknown,
  page: number,
  setPage: (page: number) => void,
): void {
  useEffect(() => {
    if (page > 0 && error instanceof ApiError && error.status === 400) setPage(0);
  }, [error, page, setPage]);
}
