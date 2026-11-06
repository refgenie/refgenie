import { useParams } from 'react-router-dom';
import { useAssetClass } from '../hooks/queries/useAssetClasses';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useUiConfig } from '../hooks/useUiConfig';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { MiniHero } from '../components/layout/MiniHero';
import { DescriptionList } from '../components/common/DescriptionList';
import { ErrorState, LoadingState } from '../components/common/states';
import { ServingModeBadges } from '../components/assets/ServingModeBadges';
import { AssetClassRecipeList } from '../components/assetClasses/AssetClassRecipeList';
import { AssetClassGenomeTable } from '../components/assetClasses/AssetClassGenomeTable';

export function AssetClassPage() {
  const { id } = useParams<{ id: string }>();
  const config = useUiConfig();
  const [state, actions] = useSearchParamsState();
  const numericId = Number(id);
  const query = useAssetClass(Number.isFinite(numericId) ? numericId : undefined);

  // No early return: the head renders in every state.
  const assetClass = query.data;

  return (
    <div className="flex flex-col gap-8">
      <MiniHero
        title={assetClass ? assetClass.name : 'Asset class'}
        documentTitle={assetClass ? undefined : false}
        breadcrumbs={
          <Breadcrumbs
            items={[
              { label: 'Asset classes', to: '/asset-classes' },
              ...(assetClass ? [{ label: assetClass.name }] : []),
            ]}
          />
        }
      />

      {query.isPending && <LoadingState label="Loading asset class" />}
      {query.error && <ErrorState error={query.error} subject={id} apiBase={config.api_base} />}

      {assetClass && (
        <>
          <DescriptionList
            items={[
              { term: 'Version', value: assetClass.version },
              { term: 'Description', value: assetClass.description ?? '' },
              {
                term: 'Serving modes',
                value: <ServingModeBadges modes={assetClass.serving_modes} />,
              },
            ]}
          />

          {/* `id` may be null on the wire; render nothing rather than filter on null. */}
          {assetClass.id !== null && (
            <>
              <AssetClassRecipeList
                assetClassId={assetClass.id}
                assetClassName={assetClass.name}
              />
              <AssetClassGenomeTable
                assetClassId={assetClass.id}
                assetClassName={assetClass.name}
                page={state.page}
                pageSize={state.pageSize}
                onPageChange={actions.setPage}
              />
            </>
          )}
        </>
      )}
    </div>
  );
}
