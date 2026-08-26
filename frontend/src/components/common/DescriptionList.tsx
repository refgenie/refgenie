import { Fragment } from 'react';
import type { ReactNode } from 'react';

export interface DescriptionItem {
  term: string;
  value: ReactNode;
}

export interface DescriptionListProps {
  items: DescriptionItem[];
  /** Drop rows whose value is null/undefined/'' instead of showing NA. */
  hideEmpty?: boolean;
}

function isEmpty(value: ReactNode): boolean {
  return value === null || value === undefined || value === '';
}

/**
 * The `NA` fallback for empty values is built in here — it is what the retired
 * `safe_show_asset_attr` Jinja macro did.
 */
export function DescriptionList({ items, hideEmpty = false }: DescriptionListProps) {
  const rows = hideEmpty ? items.filter((item) => !isEmpty(item.value)) : items;
  if (rows.length === 0) return null;
  return (
    <dl className="rg-kv">
      {rows.map((item) => (
        <Fragment key={item.term}>
          <dt className="rg-kv__term">{item.term}</dt>
          <dd className="rg-kv__value">
            {isEmpty(item.value) ? <span className="rg-muted">NA</span> : item.value}
          </dd>
        </Fragment>
      ))}
    </dl>
  );
}
