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

  // The head renders in every state, so nothing may early-return past it: a
  // loading or missing class still gets a heading, a breadcrumb out, and (via
  // MiniHero) a correct document title. No explainer here -- `/asset-classes`
  // above defines the term, and the breadcrumb links straight back to it.
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

          {/*
            Both children filter the server by asset class NAME and then guard on
            the route's numeric ID. `assetClass.id` is `number | null` on the
            wire; a null id means the record cannot be the target of this route,
            so render nothing rather than a list filtered against `null`.
          */}
          {assetClass.id !== null && (
            <>
              <AssetClassRecipeList
                assetClassId={assetClass.id}
                assetClassName={assetClass.name}
              />
              <AssetClassGenomeTable
                assetClassId={assetClass.id}
                assetClassName={assetClass.name}
                offset={state.offset}
                limit={state.limit}
                onOffsetChange={actions.setOffset}
              />
            </>
          )}
        </>
      )}
    </div>
  );
}
