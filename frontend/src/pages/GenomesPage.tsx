import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useGenomes } from '../hooks/queries/useGenomes';
import { useSummary } from '../hooks/queries/useServerInfo';
import { useUiConfig } from '../hooks/useUiConfig';
import { useCapability } from '../hooks/useCapability';
import { GenomeTable } from '../components/genomes/GenomeTable';
import { SearchBox } from '../components/common/SearchBox';
import { Pagination } from '../components/common/Pagination';
import { EmptyState } from '../components/common/states';
import { GENOME_SEARCH_FIELDS } from '../services/resources/genomes';

export function GenomesPage() {
  const config = useUiConfig();
  const [state, actions] = useSearchParamsState();
  const canReadArchives = useCapability('archives');

  const query = useGenomes({
    q: state.q || undefined,
    searchFields: state.fields.length ? state.fields : undefined,
    operator: state.operator,
    offset: state.offset,
    limit: state.limit,
  });

  // Fills the <h3>Summary</h3> heading that was empty in the Jinja page.
  const summary = useSummary({ enabled: canReadArchives });

  return (
    <div className="flex flex-col gap-6">
      <header>
        <h1 className="text-3xl font-bold mb-2">Genomes</h1>
        <p className="rg-muted">
          {config.service_name} serves reference genome assets through the refgenie API.
        </p>
      </header>

      {summary.data && (
        <dl className="grid grid-cols-1 sm:grid-cols-3 gap-4">
          {(
            [
              ['Genomes', summary.data.genomes],
              ['Asset groups', summary.data.asset_groups],
              ['Assets', summary.data.assets],
            ] as const
          ).map(([term, value]) => (
            <div className="rg-card p-4" key={term}>
              <dt className="rg-kv__term">{term}</dt>
              <dd className="text-2xl font-semibold">{value}</dd>
            </div>
          ))}
        </dl>
      )}

      <SearchBox
        value={state.q}
        onChange={actions.setQuery}
        fields={GENOME_SEARCH_FIELDS}
        selectedFields={state.fields}
        onFieldsChange={actions.setFields}
        operator={state.operator}
        onOperatorChange={actions.setOperator}
        placeholder="Search genomes…"
        label="Search genomes"
      />

      <GenomeTable
        genomes={query.data?.items}
        loading={query.isPending}
        error={query.error}
        onRetry={() => query.refetch()}
        empty={
          <EmptyState
            query={state.q || undefined}
            message="No genomes yet. Use `refgenie pull` or `refgenie build` to add one."
            onClearSearch={actions.clearSearch}
          />
        }
      />

      <Pagination pagination={query.data?.pagination} onOffsetChange={actions.setOffset} />
    </div>
  );
}
