import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useAssetClasses } from '../hooks/queries/useAssetClasses';
import { MiniHero } from '../components/layout/MiniHero';
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
      <MiniHero
        title="Asset classes"
        actions={<ActionBar slot="asset-class" />}
        lede={
          <>
            An asset class is the definition of one kind of genome file bundle — a FASTA, a
            bowtie2 index, a set of annotations. It names the files that bundle must contain, so
            anything built to that class can be used the same way whichever genome it came from.
          </>
        }
      />

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
