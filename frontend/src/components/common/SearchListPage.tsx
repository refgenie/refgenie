import type { ReactNode } from 'react';
import { useSearchParamsState } from '../../hooks/useSearchParamsState';
import { useResetPastEndPage } from '../../hooks/useResetPastEndPage';
import { MiniHero } from '../layout/MiniHero';
import { DataTable } from './DataTable';
import { SearchBox } from './SearchBox';
import { Pagination } from './Pagination';
import { EmptyState } from './states';
import type { Column } from './DataTable';
import type { ListParams, PageResponse } from '../../types/pagination';

export interface SearchListPageProps<T> {
  title: string;
  lede: ReactNode;
  actions?: ReactNode;
  useList: (p: ListParams) => {
    data: PageResponse<T> | undefined;
    isPending: boolean;
    error: unknown;
    refetch: () => unknown;
  };
  searchFields: readonly string[];
  columns: Array<Column<T>>;
  rowKey: (row: T, index: number) => string;
  emptyMessage: string;
  compact?: boolean;
  children?: ReactNode;
}

/** A catalog list page: hero, URL-backed search, table and pagination. */
export function SearchListPage<T>({
  title,
  lede,
  actions,
  useList,
  searchFields,
  columns,
  rowKey,
  emptyMessage,
  compact,
  children,
}: SearchListPageProps<T>) {
  const [state, search] = useSearchParamsState();
  const query = useList({
    q: state.q || undefined,
    searchFields: state.fields.length ? state.fields : undefined,
    operator: state.operator,
    page: state.page,
    pageSize: state.pageSize,
  });
  useResetPastEndPage(query.error, state.page, search.setPage);
  const noun = title.toLowerCase();

  return (
    <div className="flex flex-col gap-6">
      <MiniHero title={title} lede={lede} actions={actions} />
      {children}
      <SearchBox
        value={state.q}
        onChange={search.setQuery}
        fields={searchFields}
        selectedFields={state.fields}
        onFieldsChange={search.setFields}
        operator={state.operator}
        onOperatorChange={search.setOperator}
        placeholder={`Search ${noun}…`}
        label={`Search ${noun}`}
      />
      <DataTable
        caption={title}
        columns={columns}
        rows={query.data?.results}
        rowKey={rowKey}
        compact={compact}
        loading={query.isPending}
        error={query.error}
        onRetry={() => query.refetch()}
        empty={
          <EmptyState
            query={state.q || undefined}
            message={emptyMessage}
            onClearSearch={search.clearSearch}
          />
        }
      />
      <Pagination pagination={query.data?.pagination} onPageChange={search.setPage} />
    </div>
  );
}
