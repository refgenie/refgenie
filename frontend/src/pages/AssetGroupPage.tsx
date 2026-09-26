import { Link, useParams } from 'react-router-dom';
import { ScopedLink } from '../components/common/ScopedLink';
import { useAssetGroup } from '../hooks/queries/useAssetGroups';
import { useAssetClassIndex } from '../hooks/queries/useAssetClasses';
import { useAssets } from '../hooks/queries/useAssets';
import { useGenomeIndex } from '../hooks/queries/useGenomes';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
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
  // Paged, not a hardcoded `limit: 200` that silently truncates a large group.
  // `asset_group_id` is the only AssetGroup-joining filter here, which is what
  // `ListAssetsParams` requires — any second one is a server 500.
  const assets = useAssets(
    {
      asset_group_id: Number.isFinite(numericId) ? numericId : undefined,
      offset: state.offset,
      limit: state.limit,
    },
    { enabled: Number.isFinite(numericId) },
  );
  // The record carries a bare `genome_digest` and a bare `asset_class_id`;
  // both indexes exist so this page can name them instead of printing them.
  const genomes = useGenomeIndex();
  const classIndex = useAssetClassIndex();

  // The head renders in every state, so nothing may early-return past it: a
  // loading or missing group still gets a heading, a breadcrumb out, and (via
  // MiniHero) a correct document title.
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
        // The tab title qualifies the group name with its genome, because
        // `fasta` alone names a group on every genome in the instance.
        // `false` while the record loads: never title a tab "Asset group".
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
        // The second detail-page exception to "list pages define the concept":
        // there has never been an `/asset-groups` list page, so nothing above
        // this page defines "asset group".
        lede={
          <>
            An asset group holds every build of one asset class for one genome — all of this
            genome's bowtie2 indexes, say, across tool versions. One member is marked the
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
                  <Link className="rg-link" to={`/asset-classes/${groupData.asset_class_id}`}>
                    {classIndex.data?.byId.get(groupData.asset_class_id)?.name ??
                      `#${groupData.asset_class_id}`}
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
            <Pagination
              pagination={assets.data?.pagination}
              onOffsetChange={actions.setOffset}
            />
          </section>
        </>
      )}
    </div>
  );
}
