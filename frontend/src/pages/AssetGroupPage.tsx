import { useParams } from 'react-router-dom';
import { ScopedLink } from '../components/common/ScopedLink';
import { useAssetGroup } from '../hooks/queries/useAssetGroups';
import { useAssetClassIndex } from '../hooks/queries/useAssetClasses';
import { useAssets } from '../hooks/queries/useAssets';
import { useGenomeIndex } from '../hooks/queries/useGenomes';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useResetPastEndPage } from '../hooks/useResetPastEndPage';
import { useUiConfig } from '../hooks/useUiConfig';
import { genomeLabel } from '../utils/genomes';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { MiniHero } from '../components/layout/MiniHero';
import { DescriptionList } from '../components/common/DescriptionList';
import { DigestChip } from '../components/common/DigestChip';
import { ErrorState, LoadingState } from '../components/common/states';
import { Pagination } from '../components/common/Pagination';
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
  const [state, actions] = useSearchParamsState();
  const numericId = Number(id);
  const group = useAssetGroup(Number.isFinite(numericId) ? numericId : undefined);
  // `ListAssetsParams` allows only one AssetGroup-joining filter; a second one is a server 500.
  const assets = useAssets(
    {
      asset_group_id: Number.isFinite(numericId) ? numericId : undefined,
      page: state.page,
      pageSize: state.pageSize,
    },
    { enabled: Number.isFinite(numericId) },
  );
  useResetPastEndPage(assets.error, state.page, actions.setPage);
  // Indexes to show names for the record's genome digest and asset class id.
  const genomes = useGenomeIndex();
  const classIndex = useAssetClassIndex();

  // No early return: the head renders in every state.
  const groupData = group.data;
  const curationColumns: Array<Column<AssetResponse>> = groupData
    ? [
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
      ]
    : [];

  return (
    <div className="flex flex-col gap-8">
      <MiniHero
        title={groupData ? groupData.name : 'Asset group'}
        // Qualify with the genome (`fasta` exists on every genome); `false` while loading.
        documentTitle={
          groupData
            ? `${groupData.name} · ${genomeLabel(genomes.data, groupData.genome_digest)}`
            : false
        }
        breadcrumbs={
          <Breadcrumbs
            items={[
              { label: 'Genomes', to: '/genomes' },
              ...(groupData
                ? [
                    {
                      label: genomeLabel(genomes.data, groupData.genome_digest),
                      to: `/genomes/${groupData.genome_digest}`,
                    },
                    { label: groupData.name },
                  ]
                : []),
            ]}
          />
        }
        // No `/asset-groups` list page exists, so the lede defines the term here.
        lede={
          <>
            An asset group holds every build of one asset class for one genome: all of this
            genome's bowtie2 indexes across tool versions, for example. One member is marked the
            default, and that is what a bare{' '}
            <code className="rg-code rg-code--inline">{'<genome>/<group>'}</code> registry path
            resolves to.
          </>
        }
      />

      {group.isPending && <LoadingState label="Loading asset group" />}
      {group.error && <ErrorState error={group.error} subject={id} apiBase={config.api_base} />}

      {groupData && (
        <>
          <DescriptionList
            items={[
              { term: 'Description', value: groupData.description ?? '' },
              {
                term: 'Genome',
                value: (
                  <span className="flex items-center gap-2 flex-wrap">
                    <ScopedLink className="rg-link" to={`/genomes/${groupData.genome_digest}`}>
                      {genomeLabel(genomes.data, groupData.genome_digest)}
                    </ScopedLink>
                    <DigestChip digest={groupData.genome_digest} />
                  </span>
                ),
              },
              {
                term: 'Asset class',
                value: (
                  <ScopedLink className="rg-link" to={`/asset-classes/${groupData.asset_class_id}`}>
                    {classIndex.data?.byId.get(groupData.asset_class_id)?.name ??
                      `#${groupData.asset_class_id}`}
                  </ScopedLink>
                ),
              },
            ]}
          />

          <section>
            <h2 className="text-xl font-semibold mb-4">Assets</h2>
            <AssetTable
              assets={assets.data?.results}
              loading={assets.isPending}
              error={assets.error}
              onRetry={() => assets.refetch()}
              extraColumns={curationColumns}
            />
            <Pagination
              pagination={assets.data?.pagination}
              onPageChange={actions.setPage}
            />
          </section>
        </>
      )}
    </div>
  );
}
