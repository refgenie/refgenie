import { Link, useParams } from 'react-router-dom';
import { useAssetGroup } from '../hooks/queries/useAssetGroups';
import { useAssets } from '../hooks/queries/useAssets';
import { useUiConfig } from '../hooks/useUiConfig';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { DescriptionList } from '../components/common/DescriptionList';
import { ErrorState, LoadingState } from '../components/common/states';
import { AssetTable } from '../components/assets/AssetTable';
import { SetDefaultControl } from '../components/actions/SetDefaultControl';
import { DeleteAssetButton } from '../components/actions/DeleteAssetButton';
import type { Column } from '../components/common/DataTable';
import type { AssetResponse } from '../types/api';

function isDefaultAsset(asset: AssetResponse): boolean {
  return asset.is_default === true;
}

export function AssetGroupPage() {
  const { id } = useParams<{ id: string }>();
  const config = useUiConfig();
  const numericId = Number(id);
  const group = useAssetGroup(Number.isFinite(numericId) ? numericId : undefined);
  const assets = useAssets(
    { asset_group_id: Number.isFinite(numericId) ? numericId : undefined, limit: 200 },
    { enabled: Number.isFinite(numericId) },
  );

  if (group.isPending) return <LoadingState label="Loading asset group" />;
  if (group.error) {
    return <ErrorState error={group.error} subject={id} apiBase={config.api_base} />;
  }
  if (!group.data) return null;

  const groupData = group.data;
  const curationColumns: Array<Column<AssetResponse>> = [
    {
      key: 'default',
      header: 'Default',
      render: (asset) => (
        <SetDefaultControl
          genomeDigest={groupData.genome_digest}
          assetGroupName={groupData.name}
          assetName={asset.name}
          isDefault={isDefaultAsset(asset)}
        />
      ),
    },
    {
      key: 'actions',
      header: '',
      render: (asset) =>
        asset.digest ? (
          <DeleteAssetButton
            compact
            digest={asset.digest}
            registryPath={`${groupData.name}:${asset.name}`}
            seekKeys={asset.seek_keys}
          />
        ) : null,
    },
  ];

  return (
    <div className="flex flex-col gap-8">
      <Breadcrumbs
        items={[
          { label: 'Genomes', to: '/genomes' },
          { label: 'Genome', to: `/genomes/${group.data.genome_digest}` },
          { label: group.data.name },
        ]}
      />
      <h1 className="text-3xl font-bold">{group.data.name}</h1>

      <DescriptionList
        items={[
          { term: 'Description', value: group.data.description ?? '' },
          {
            term: 'Genome',
            value: (
              <Link className="rg-link" to={`/genomes/${group.data.genome_digest}`}>
                {group.data.genome_digest}
              </Link>
            ),
          },
          {
            term: 'Asset class',
            value: (
              <Link className="rg-link" to={`/asset-classes/${group.data.asset_class_id}`}>
                #{group.data.asset_class_id}
              </Link>
            ),
          },
        ]}
      />

      <section>
        <h2 className="text-xl font-semibold mb-4">Assets</h2>
        <AssetTable
          assets={assets.data?.items}
          loading={assets.isPending}
          error={assets.error}
          onRetry={() => assets.refetch()}
          extraColumns={curationColumns}
        />
      </section>
    </div>
  );
}
