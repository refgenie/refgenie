/**
 * Search and pagination state kept in the URL, so filtered views are shareable
 * and survive reload. Filter changes replace the history entry (SearchBox
 * debounces keystrokes). The URL `page` is 1-based; `state.page` is 0-based for the API.
 */

import { useCallback, useMemo } from 'react';
import { useSearchParams } from 'react-router-dom';
import { DEFAULT_PAGE_SIZE } from '../types/pagination';
import type { SearchOperator } from '../types/pagination';

export interface ListViewState {
  q: string;
  fields: string[];
  operator: SearchOperator;
  /** 0-based, ready for the API. */
  page: number;
  pageSize: number;
}

export interface ListViewActions {
  setQuery: (q: string) => void;
  setFields: (fields: string[]) => void;
  setOperator: (operator: SearchOperator) => void;
  /** Takes the 0-based API page. */
  setPage: (page: number) => void;
  clearSearch: () => void;
}

const OPERATORS = new Set<SearchOperator>(['eq', 'contains', 'starts_with', 'ends_with']);

function parseOperator(raw: string | null): SearchOperator {
  return raw && OPERATORS.has(raw as SearchOperator) ? (raw as SearchOperator) : 'contains';
}

function parseInteger(raw: string | null, fallback: number): number {
  // Number(null) === 0, so we must explicitly check for null/empty before conversion
  if (!raw) return fallback;
  const parsed = Number(raw);
  return Number.isFinite(parsed) && parsed >= 0 ? Math.floor(parsed) : fallback;
}

export function useSearchParamsState(): [ListViewState, ListViewActions] {
  const [params, setParams] = useSearchParams();

  const state = useMemo<ListViewState>(() => {
    const rawFields = params.get('fields');
    return {
      q: params.get('q') ?? '',
      fields: rawFields ? rawFields.split(',').filter(Boolean) : [],
      operator: parseOperator(params.get('op')),
      page: Math.max(0, parseInteger(params.get('page'), 1) - 1),
      pageSize: Math.max(1, parseInteger(params.get('page_size'), DEFAULT_PAGE_SIZE)),
    };
  }, [params]);

  const update = useCallback(
    (patch: Record<string, string | undefined>, resetPage: boolean) => {
      setParams(
        (previous) => {
          const next = new URLSearchParams(previous);
          for (const [key, value] of Object.entries(patch)) {
            if (value === undefined || value === '') next.delete(key);
            else next.set(key, value);
          }
          if (resetPage) next.delete('page');
          return next;
        },
        { replace: true },
      );
    },
    [setParams],
  );

  const actions = useMemo<ListViewActions>(
    () => ({
      // Any query change resets to the first page.
      setQuery: (q) => update({ q }, true),
      setFields: (fields) => update({ fields: fields.join(',') }, true),
      setOperator: (operator) => update({ op: operator }, true),
      // Paging is the only change that keeps the page.
      setPage: (page) => update({ page: page > 0 ? String(page + 1) : undefined }, false),
      clearSearch: () => update({ q: undefined, fields: undefined, op: undefined }, true),
    }),
    [update],
  );

  return [state, actions];
}
