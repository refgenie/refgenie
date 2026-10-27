import { useState } from 'react';
import { Link } from 'react-router-dom';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useAliases } from '../hooks/queries/useAliases';
import { useCapability } from '../hooks/useCapability';
import { AddAliasModal } from '../components/actions/AddAliasModal';
import { MiniHero } from '../components/layout/MiniHero';
import { DataTable } from '../components/common/DataTable';
import { DigestChip } from '../components/common/DigestChip';
import { SearchBox } from '../components/common/SearchBox';
import { Pagination } from '../components/common/Pagination';
import { EmptyState } from '../components/common/states';
import { ALIAS_SEARCH_FIELDS } from '../services/resources/aliases';
import type { Column } from '../components/common/DataTable';
import type { AliasPublic } from '../types/api';

export function AliasesPage() {
  const canWrite = useCapability('aliases_write');
  const [addOpen, setAddOpen] = useState(false);
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
      <MiniHero
        title="Aliases"
        lede={
          <>
            An alias is a human-readable name — <code className="rg-code rg-code--inline">hg38</code>,{' '}
            <code className="rg-code rg-code--inline">mm10</code> — that resolves to a genome
            digest. You type the alias; refgenie stores and compares everything by the digest,
            and one genome can answer to several aliases.
          </>
        }
        actions={
          canWrite && (
            <button
              type="button"
              className="rg-btn rg-btn--primary"
              onClick={() => setAddOpen(true)}
            >
              Add alias
            </button>
          )
        }
      />

      {canWrite && <AddAliasModal isOpen={addOpen} onClose={() => setAddOpen(false)} />}

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
