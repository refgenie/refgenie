import type { ReactNode } from 'react';
import { cn } from '../../utils/cn';
import { EmptyState, ErrorState, LoadingState } from './states';

export interface Column<T> {
  key: string;
  header: string;
  render: (row: T) => ReactNode;
  align?: 'left' | 'right';
  width?: string;
}

export interface DataTableProps<T> {
  /** Screen-reader caption; every table has one. */
  caption: string;
  columns: Array<Column<T>>;
  rows: T[] | undefined;
  rowKey: (row: T, index: number) => string;
  loading?: boolean;
  error?: unknown;
  /** Rendered when `rows` is empty and there is no error. */
  empty?: ReactNode;
  /**
   * Optional whole-row click. The primary cell must ALSO contain a real
   * <a>/<Link> so keyboard users and middle-click work.
   */
  onRowClick?: (row: T) => void;
  onRetry?: () => void;
}

export function DataTable<T>({
  caption,
  columns,
  rows,
  rowKey,
  loading,
  error,
  empty,
  onRowClick,
  onRetry,
}: DataTableProps<T>) {
  if (loading) return <LoadingState rows={5} label={`Loading ${caption}`} />;
  if (error) return <ErrorState error={error} onRetry={onRetry} />;
  if (!rows || rows.length === 0) {
    return <>{empty ?? <EmptyState message={`No ${caption.toLowerCase()}.`} />}</>;
  }

  return (
    <div className="rg-table__wrap">
      <table className="rg-table">
        <caption className="sr-only">{caption}</caption>
        <thead className="rg-table__head">
          <tr>
            {columns.map((column) => (
              <th
                key={column.key}
                scope="col"
                className={cn(
                  'rg-table__cell',
                  column.align === 'right' && 'rg-table__cell--numeric',
                )}
              >
                {column.header}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map((row, index) => (
            <tr
              key={rowKey(row, index)}
              className={cn('rg-table__row', onRowClick && 'rg-table__row--clickable')}
              onClick={onRowClick ? () => onRowClick(row) : undefined}
            >
              {columns.map((column) => (
                <td
                  key={column.key}
                  className={cn(
                    'rg-table__cell',
                    column.align === 'right' && 'rg-table__cell--numeric',
                  )}
                >
                  {column.render(row)}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
