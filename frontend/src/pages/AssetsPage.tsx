import { useSearchParams } from 'react-router-dom';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useAssets } from '../hooks/queries/useAssets';
import { AssetTable } from '../components/assets/AssetTable';
import { SearchBox } from '../components/common/SearchBox';
import { Pagination } from '../components/common/Pagination';
import { EmptyState } from '../components/common/states';
import { ASSET_SEARCH_FIELDS } from '../services/resources/assets';

export function AssetsPage() {
  const [state, actions] = useSearchParamsState();
  const [params] = useSearchParams();
  // Deep links from the genome table carry a genome filter.
  const genomeDigest = params.get('genome_digest') ?? undefined;

  const query = useAssets({
    genome_digest: genomeDigest,
    q: state.q || undefined,
    searchFields: state.fields.length ? state.fields : undefined,
    operator: state.operator,
    offset: state.offset,
    limit: state.limit,
  });

  return (
    <div className="flex flex-col gap-6">
      <header>
        <h1 className="text-3xl font-bold mb-2">Assets</h1>
        {genomeDigest && (
          <p className="rg-muted text-sm">
            Filtered to genome <code className="rg-code rg-code--inline">{genomeDigest}</code>
          </p>
        )}
      </header>

      <SearchBox
        value={state.q}
        onChange={actions.setQuery}
        fields={ASSET_SEARCH_FIELDS}
        selectedFields={state.fields}
        onFieldsChange={actions.setFields}
        operator={state.operator}
        onOperatorChange={actions.setOperator}
        placeholder="Search assets…"
        label="Search assets"
      />

      <AssetTable
        assets={query.data?.items}
        loading={query.isPending}
        error={query.error}
        onRetry={() => query.refetch()}
        empty={
          <EmptyState
            query={state.q || undefined}
            message="No assets yet. Use `refgenie pull` or `refgenie build` to add one."
            onClearSearch={actions.clearSearch}
          />
        }
      />

      <Pagination pagination={query.data?.pagination} onOffsetChange={actions.setOffset} />
    </div>
  );
}
