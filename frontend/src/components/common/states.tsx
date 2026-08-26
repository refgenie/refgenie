/**
 * The three mandatory states for every data-driven view. Never a blank screen,
 * never a silently empty list.
 */

import { ApiError } from '../../services/http';
import { CopyButton } from './CopyButton';

export interface LoadingStateProps {
  /** Render skeleton rows (tables) instead of a spinner block (panels). */
  rows?: number;
  label?: string;
}

export function LoadingState({ rows = 0, label = 'Loading' }: LoadingStateProps) {
  if (rows > 0) {
    return (
      <div className="p-4 flex flex-col gap-3" aria-busy="true">
        <span className="sr-only">{label}</span>
        {Array.from({ length: rows }).map((_, i) => (
          <span className="rg-skeleton" key={i} />
        ))}
      </div>
    );
  }
  return (
    <div className="rg-state rg-state--loading" aria-busy="true">
      {label}…
    </div>
  );
}

export interface ErrorStateProps {
  error: unknown;
  /** Shown in the 404 message and given a copy button. */
  subject?: string;
  /** API base, named in the network-failure message. */
  apiBase?: string;
  onRetry?: () => void;
}

export function ErrorState({ error, subject, apiBase, onRetry }: ErrorStateProps) {
  const api = error instanceof ApiError ? error : undefined;
  const detail = api?.detail ?? (error instanceof Error ? error.message : String(error));

  let body: React.ReactNode = detail;
  if (api?.isNetworkError) {
    body = `Cannot reach the refgenie API at ${apiBase ?? api.url}. Is the server running?`;
  } else if (api?.isNotFound) {
    body = (
      <span className="flex items-center gap-2 flex-wrap">
        <span>Not found{subject ? ':' : '.'}</span>
        {subject && <code className="rg-code rg-code--inline">{subject}</code>}
        {subject && <CopyButton value={subject} label="Copy identifier" />}
      </span>
    );
  }

  return (
    <div className="rg-state rg-state--error" role="alert">
      <p className="font-semibold mb-1">Request failed{api ? ` (${api.status})` : ''}</p>
      <div className="text-sm">{body}</div>
      {onRetry && api && api.status >= 500 && (
        <button type="button" className="rg-btn rg-btn--sm mt-4" onClick={onRetry}>
          Retry
        </button>
      )}
    </div>
  );
}

export interface EmptyStateProps {
  /** The active search term, if any. Distinguishes "no results" from "no data". */
  query?: string;
  /** Shown when there is no data at all. Text only: buttons are the manage UI. */
  message: string;
  onClearSearch?: () => void;
}

export function EmptyState({ query, message, onClearSearch }: EmptyStateProps) {
  if (query) {
    return (
      <div className="rg-state rg-state--empty">
        <p className="mb-4">No results match “{query}”.</p>
        {onClearSearch && (
          <button type="button" className="rg-btn rg-btn--sm" onClick={onClearSearch}>
            Clear search
          </button>
        )}
      </div>
    );
  }
  return <div className="rg-state rg-state--empty">{message}</div>;
}
