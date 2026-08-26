import { Link } from 'react-router-dom';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useAliases } from '../hooks/queries/useAliases';
import { DataTable } from '../components/common/DataTable';
import { DigestChip } from '../components/common/DigestChip';
import { SearchBox } from '../components/common/SearchBox';
import { Pagination } from '../components/common/Pagination';
import { EmptyState } from '../components/common/states';
import { ALIAS_SEARCH_FIELDS } from '../services/resources/aliases';
import type { Column } from '../components/common/DataTable';
import type { AliasPublic } from '../types/api';

export function AliasesPage() {
  const [state, actions] = useSearchParamsState();
  const query = useAliases({
    q: state.q || undefined,
    searchFields: state.fields.length ? state.fields : undefined,
    operator: state.operator,
    offset: state.offset,
    limit: state.limit,
  });

  const columns: Array<Column<AliasPublic>> = [
    {
      key: 'name',
      header: 'Alias',
      render: (alias) => (
        <Link className="rg-link font-medium" to={`/genomes/${alias.genome_digest}`}>
          {alias.name}
        </Link>
      ),
    },
    {
      key: 'digest',
      header: 'Genome digest',
      render: (alias) => <DigestChip digest={alias.genome_digest} />,
    },
  ];

  return (
    <div className="flex flex-col gap-6">
      <h1 className="text-3xl font-bold">Aliases</h1>

      <SearchBox
        value={state.q}
        onChange={actions.setQuery}
        fields={ALIAS_SEARCH_FIELDS}
        selectedFields={state.fields}
        onFieldsChange={actions.setFields}
        operator={state.operator}
        onOperatorChange={actions.setOperator}
        placeholder="Search aliases…"
        label="Search aliases"
      />

      <DataTable
        caption="Aliases"
        columns={columns}
        rows={query.data?.items}
        rowKey={(alias) => alias.name}
        loading={query.isPending}
        error={query.error}
        onRetry={() => query.refetch()}
        empty={
          <EmptyState
            query={state.q || undefined}
            message="No aliases defined."
            onClearSearch={actions.clearSearch}
          />
        }
      />

      <Pagination pagination={query.data?.pagination} onOffsetChange={actions.setOffset} />
    </div>
  );
}
