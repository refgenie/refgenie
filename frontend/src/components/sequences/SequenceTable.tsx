// Ported from repos/refget/frontend/src/components/SequenceTable.tsx @ c0d3e09; keep the two in step.

import { useMemo, useState } from 'react';
import type { ReactNode } from 'react';
import { Link } from 'react-router-dom';
import { DigestChip } from '../common/DigestChip';
import { EmptyState } from '../common/states';
import { Pagination } from '../common/Pagination';
import { SearchBox } from '../common/SearchBox';
import { cn } from '../../utils/cn';
import { preferredAlias } from './aliases';
import type { SequenceRow } from './types';
import { isExtractable } from '../../services/regionExtract';

export const SEQUENCE_PAGE_SIZE = 50;

type SortColumn = 'name' | 'length';
type SortState = { column: SortColumn; ascending: boolean } | null;

export interface SequenceTableProps {
  sequences: SequenceRow[];
  /** Shows a per-row "Extract…" link, to this href, for rows whose bases can be read. */
  extractHref?: (row: SequenceRow) => string;
  /** Shown to the right of the filter box, e.g. the alphabet summary. */
  aside?: ReactNode;
}

function matches(row: SequenceRow, needle: string): boolean {
  if (row.name.toLowerCase().includes(needle)) return true;
  return row.aliases.some((a) => a.alias.toLowerCase().includes(needle));
}

function compareRows(a: SequenceRow, b: SequenceRow, column: SortColumn): number {
  return column === 'length' ? a.length - b.length : a.name.localeCompare(b.name);
}

const NO_ALIAS =
  'No alias: no public record, such as RefSeq, has exactly these bases. Masked or edited sequences often have none.';

function AliasCell({ row }: { row: SequenceRow }) {
  const { primary, all } = preferredAlias(row.aliases);
  if (!primary) {
    // Not an error: masked or edited sequences match no public record.
    return (
      <span className="rg-muted" title={NO_ALIAS}>
        <span aria-hidden="true">–</span>
        <span className="sr-only">{NO_ALIAS}</span>
      </span>
    );
  }
  const tooltip = all.join('\n');
  return (
    <span title={tooltip}>
      {primary}
      {all.length > 1 && <span className="rg-muted"> and {all.length - 1} more</span>}
    </span>
  );
}

/**
 * A genome's sequences: filter by name or alias, sort by name or length,
 * 50 rows a page. Presentational: rows in, table out.
 */
export function SequenceTable({ sequences, extractHref, aside }: SequenceTableProps) {
  const [query, setQuery] = useState('');
  const [sort, setSort] = useState<SortState>(null);
  const [page, setPage] = useState(0);

  const rows = useMemo(() => {
    const needle = query.trim().toLowerCase();
    const filtered = needle ? sequences.filter((row) => matches(row, needle)) : sequences;
    if (!sort) return filtered;
    const sign = sort.ascending ? 1 : -1;
    return [...filtered].sort((a, b) => sign * compareRows(a, b, sort.column));
  }, [sequences, query, sort]);

  const paged = rows.slice(page * SEQUENCE_PAGE_SIZE, (page + 1) * SEQUENCE_PAGE_SIZE);

  const onQuery = (q: string) => {
    setQuery(q);
    setPage(0);
  };

  const toggleSort = (column: SortColumn) => {
    setSort((current) =>
      current?.column === column
        ? { column, ascending: !current.ascending }
        : { column, ascending: true },
    );
    setPage(0);
  };

  const ariaSort = (column: SortColumn) =>
    sort?.column === column ? (sort.ascending ? 'ascending' : 'descending') : undefined;

  const sortHeader = (column: SortColumn, label: string, numeric = false) => (
    <th
      scope="col"
      aria-sort={ariaSort(column)}
      className={cn('rg-table__cell', numeric && 'rg-table__cell--numeric')}
    >
      <button type="button" className="rg-seq-sort" onClick={() => toggleSort(column)}>
        {label}
        {sort?.column === column && <span aria-hidden="true">{sort.ascending ? '▲' : '▼'}</span>}
      </button>
    </th>
  );

  return (
    <div className="flex flex-col gap-4">
      <div className="rg-seq-toolbar">
        <SearchBox
          value={query}
          onChange={onQuery}
          label="Filter sequences"
          placeholder="Filter by name or alias…"
        />
        {aside}
      </div>
      {rows.length === 0 ? (
        <EmptyState query={query} message="No sequences." onClearSearch={() => onQuery('')} />
      ) : (
        <div className="rg-table__wrap">
          <table className="rg-table rg-table--compact rg-table--stack">
            <caption className="sr-only">Sequences</caption>
            <thead className="rg-table__head">
              <tr>
                {sortHeader('name', 'Seq name')}
                {sortHeader('length', 'Length', true)}
                <th scope="col" className="rg-table__cell">
                  Aliases
                </th>
                <th scope="col" className="rg-table__cell">
                  SHA-512/24u
                </th>
                <th scope="col" className="rg-table__cell">
                  MD5
                </th>
                <th scope="col" className="rg-table__cell">
                  Alphabet
                </th>
                {extractHref && (
                  <th scope="col" className="rg-table__cell">
                    <span className="sr-only">Actions</span>
                  </th>
                )}
              </tr>
            </thead>
            <tbody>
              {paged.map((row, index) => (
                <tr key={`${row.sha512t24u}-${row.name}-${index}`} className="rg-table__row">
                  <td className="rg-table__cell">{row.name}</td>
                  <td className="rg-table__cell rg-table__cell--numeric" data-label="Length">
                    {row.length.toLocaleString()}
                  </td>
                  <td className="rg-table__cell" data-label="Aliases">
                    <AliasCell row={row} />
                  </td>
                  <td className="rg-table__cell" data-label="SHA-512/24u">
                    <DigestChip digest={row.sha512t24u} />
                  </td>
                  <td className="rg-table__cell" data-label="MD5">
                    <DigestChip digest={row.md5} />
                  </td>
                  <td className="rg-table__cell" data-label="Alphabet">
                    {row.alphabet || (
                      <span className="rg-muted" title="Not recorded">
                        –
                      </span>
                    )}
                  </td>
                  {extractHref && (
                    <td className="rg-table__cell" data-label="">
                      {isExtractable(row.alphabet) && (
                        <Link
                          className="rg-link text-sm"
                          aria-label={`Extract a region of ${row.name}`}
                          to={extractHref(row)}
                        >
                          Extract…
                        </Link>
                      )}
                    </td>
                  )}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
      <Pagination
        pagination={{ page, page_size: SEQUENCE_PAGE_SIZE, total: rows.length }}
        onPageChange={setPage}
      />
    </div>
  );
}
