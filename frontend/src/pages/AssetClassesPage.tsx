import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useAssetClasses } from '../hooks/queries/useAssetClasses';
import { AssetClassTable } from '../components/assetClasses/AssetClassTable';
import { SearchBox } from '../components/common/SearchBox';
import { Pagination } from '../components/common/Pagination';
import { ActionBar } from '../components/common/ActionBar';
import { EmptyState } from '../components/common/states';
import { ASSET_CLASS_SEARCH_FIELDS } from '../services/resources/assetClasses';

export function AssetClassesPage() {
  const [state, actions] = useSearchParamsState();
  const query = useAssetClasses({
    q: state.q || undefined,
    searchFields: state.fields.length ? state.fields : undefined,
    operator: state.operator,
    offset: state.offset,
    limit: state.limit,
  });

  return (
    <div className="flex flex-col gap-6">
      <header className="flex items-start justify-between gap-4 flex-wrap">
        <h1 className="text-3xl font-bold">Asset classes</h1>
        <ActionBar slot="asset-class" />
      </header>

      <SearchBox
        value={state.q}
        onChange={actions.setQuery}
        fields={ASSET_CLASS_SEARCH_FIELDS}
        selectedFields={state.fields}
        onFieldsChange={actions.setFields}
        operator={state.operator}
        onOperatorChange={actions.setOperator}
        placeholder="Search asset classes…"
        label="Search asset classes"
      />

      <AssetClassTable
        assetClasses={query.data?.items}
        loading={query.isPending}
        error={query.error}
        onRetry={() => query.refetch()}
        empty={
          <EmptyState
            query={state.q || undefined}
            message="No asset classes registered."
            onClearSearch={actions.clearSearch}
          />
        }
      />

      <Pagination pagination={query.data?.pagination} onOffsetChange={actions.setOffset} />
    </div>
  );
}
