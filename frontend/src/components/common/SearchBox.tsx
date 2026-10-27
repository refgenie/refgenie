import { useEffect, useId, useState } from 'react';
import { SEARCH_OPERATORS } from '../../types/pagination';
import type { SearchOperator } from '../../types/pagination';

export interface SearchBoxProps {
  value: string;
  onChange: (q: string) => void;
  /** Server-enforced allowlist for this endpoint; a bad value is a 422. */
  fields?: readonly string[];
  selectedFields?: string[];
  onFieldsChange?: (fields: string[]) => void;
  operator?: SearchOperator;
  onOperatorChange?: (operator: SearchOperator) => void;
  placeholder?: string;
  label?: string;
}

const DEBOUNCE_MS = 300;

/**
 * Debounced search input. Writes to the URL (via the caller), never filters an
 * in-memory array: search is always server-side.
 */
export function SearchBox({
  value,
  onChange,
  fields,
  selectedFields = [],
  onFieldsChange,
  operator = 'contains',
  onOperatorChange,
  placeholder = 'Search…',
  label = 'Search',
}: SearchBoxProps) {
  const inputId = useId();
  const fieldId = useId();
  const operatorId = useId();
  const [draft, setDraft] = useState(value);

  // Re-sync when the URL changes underneath us (back button, clear-search).
  useEffect(() => setDraft(value), [value]);

  useEffect(() => {
    if (draft === value) return;
    const timer = setTimeout(() => onChange(draft), DEBOUNCE_MS);
    return () => clearTimeout(timer);
  }, [draft, value, onChange]);

  return (
    <div className="rg-search">
      <label className="sr-only" htmlFor={inputId}>
        {label}
      </label>
      <input
        id={inputId}
        className="rg-search__input"
        type="search"
        value={draft}
        placeholder={placeholder}
        onChange={(event) => setDraft(event.target.value)}
      />

      {fields && fields.length > 0 && onFieldsChange && (
        <>
          <label className="sr-only" htmlFor={fieldId}>
            Field to search
          </label>
          <select
            id={fieldId}
            className="rg-search__field-select"
            value={selectedFields[0] ?? ''}
            onChange={(event) =>
              onFieldsChange(event.target.value ? [event.target.value] : [])
            }
          >
            <option value="">All fields</option>
            {fields.map((field) => (
              <option key={field} value={field}>
                {field}
              </option>
            ))}
          </select>
        </>
      )}

      {onOperatorChange && (
        <>
          <label className="sr-only" htmlFor={operatorId}>
            Match type
          </label>
          <select
            id={operatorId}
            className="rg-search__operator"
            value={operator}
            onChange={(event) => onOperatorChange(event.target.value as SearchOperator)}
          >
            {SEARCH_OPERATORS.map((op) => (
              <option key={op} value={op}>
                {op}
              </option>
            ))}
          </select>
        </>
      )}
    </div>
  );
}
