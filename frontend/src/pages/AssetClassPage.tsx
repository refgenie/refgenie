import { useParams } from 'react-router-dom';
import { useAssetClass } from '../hooks/queries/useAssetClasses';
import { useUiConfig } from '../hooks/useUiConfig';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { DescriptionList } from '../components/common/DescriptionList';
import { ErrorState, LoadingState } from '../components/common/states';
import { ServingModeBadges } from '../components/assets/ServingModeBadges';

export function AssetClassPage() {
  const { id } = useParams<{ id: string }>();
  const config = useUiConfig();
  const numericId = Number(id);
  const query = useAssetClass(Number.isFinite(numericId) ? numericId : undefined);

  if (query.isPending) return <LoadingState label="Loading asset class" />;
  if (query.error) {
    return <ErrorState error={query.error} subject={id} apiBase={config.api_base} />;
  }
  if (!query.data) return null;

  return (
    <div className="flex flex-col gap-8">
      <Breadcrumbs
        items={[
          { label: 'Asset classes', to: '/asset-classes' },
          { label: query.data.name },
        ]}
      />
      <h1 className="text-3xl font-bold">{query.data.name}</h1>
      <DescriptionList
        items={[
          { term: 'Version', value: query.data.version },
          { term: 'Description', value: query.data.description ?? '' },
          {
            term: 'Serving modes',
            value: <ServingModeBadges modes={query.data.serving_modes} />,
          },
        ]}
      />
    </div>
  );
}
