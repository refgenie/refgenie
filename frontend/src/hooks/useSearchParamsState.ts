/**
 * Search and pagination state lives in the URL, not in component state.
 *
 * Consequences that must hold: a filtered view is a shareable link, and a
 * reload restores the view.
 *
 * Filter changes REPLACE the history entry rather than pushing one, because
 * `SearchBox` debounces at 300 ms and pushing per debounced keystroke would
 * bury the previous page under a stack of near-identical URLs. The cost is
 * that the back button leaves the list view instead of stepping back through
 * filters; that is the intended trade, not an oversight.
 *
 * This hook never filters an in-memory array. Pages generally do not either —
 * `GenomesPage` is the one exception, because no endpoint has a `has_assets`
 * parameter and 675 of 701 genomes have nothing to show.
 */

import { useCallback, useMemo } from 'react';
import { useSearchParams } from 'react-router-dom';
import { DEFAULT_PAGE_SIZE } from '../types/pagination';
import type { SearchOperator } from '../types/pagination';

export interface ListViewState {
  q: string;
  fields: string[];
  operator: SearchOperator;
  offset: number;
  limit: number;
}

export interface ListViewActions {
  setQuery: (q: string) => void;
  setFields: (fields: string[]) => void;
  setOperator: (operator: SearchOperator) => void;
  setOffset: (offset: number) => void;
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
      offset: parseInteger(params.get('offset'), 0),
      limit: Math.max(1, parseInteger(params.get('limit'), DEFAULT_PAGE_SIZE)),
    };
  }, [params]);

  const update = useCallback(
    (patch: Record<string, string | undefined>, resetOffset: boolean) => {
      setParams(
        (previous) => {
          const next = new URLSearchParams(previous);
          for (const [key, value] of Object.entries(patch)) {
            if (value === undefined || value === '') next.delete(key);
            else next.set(key, value);
          }
          if (resetOffset) next.delete('offset');
          return next;
        },
        { replace: true },
      );
    },
    [setParams],
  );

  const actions = useMemo<ListViewActions>(
    () => ({
      // Any query change resets pagination: page 4 of the old result set is
      // meaningless against the new one.
      setQuery: (q) => update({ q }, true),
      setFields: (fields) => update({ fields: fields.join(',') }, true),
      setOperator: (operator) => update({ op: operator }, true),
      // Paging is the one change that keeps the offset. If a `setLimit` is ever
      // added it must pass `true`: page 7 at 10 rows is row 61, but page 7 at
      // 100 rows is row 601, and the user did not ask to go there.
      setOffset: (offset) => update({ offset: offset > 0 ? String(offset) : undefined }, false),
      clearSearch: () => update({ q: undefined, fields: undefined, op: undefined }, true),
    }),
    [update],
  );

  return [state, actions];
}
