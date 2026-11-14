import { DataTable } from '../common/DataTable';
import { ScopedLink } from '../common/ScopedLink';
import { DigestChip } from '../common/DigestChip';
import { DownloadLink } from '../common/DownloadLink';
import { FileSize } from '../common/FileSize';
import { DeleteAssetButton } from '../actions/DeleteAssetButton';
import { BridgePullButton } from '../bridge/BridgePullButton';
import { ServingModeBadges } from '../assets/ServingModeBadges';
import { archiveDownloadUrl } from '../../services/resources/serverInfo';
import { useApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import type { Column } from '../common/DataTable';
import type { ArchiveRecord, AssetResponse } from '../../types/api';

export interface GenomeAssetTableProps {
  assets: AssetResponse[] | undefined;
  archives: ArchiveRecord[] | undefined;
  loading?: boolean;
  error?: unknown;
  empty?: React.ReactNode;
  onRetry?: () => void;
}

interface JoinedAsset {
  asset: AssetResponse;
  archive: ArchiveRecord | undefined;
}

/**
 * The genome's assets, grouped by asset group. Built from `/assets` and
 * left-joined to archives, so file-mode assets with no tarball still get a row.
 */
export function GenomeAssetTable({
  assets,
  archives,
  loading,
  error,
  empty,
  onRetry,
}: GenomeAssetTableProps) {
  const client = useApiClient();
  const canDownloadArchives = useCapability('archives');

  const rows: JoinedAsset[] | undefined = assets
    ?.map((asset) => ({
      asset,
      archive: archives?.find((archive) => archive.asset_digest === asset.digest),
    }))
    .sort((a, b) => {
      const group = (a.asset.asset_group_name ?? '').localeCompare(
        b.asset.asset_group_name ?? '',
      );
      return group !== 0 ? group : a.asset.name.localeCompare(b.asset.name);
    });

  const columns: Array<Column<JoinedAsset>> = [
    {
      key: 'name',
      header: 'Asset',
      render: ({ asset }) =>
        asset.digest ? (
          <ScopedLink className="rg-link font-medium" to={`/assets/${asset.digest}`}>
            {asset.asset_group_name ? `${asset.asset_group_name}:` : ''}
            {asset.name}
          </ScopedLink>
        ) : (
          <span>{asset.name}</span>
        ),
    },
    {
      key: 'description',
      // The asset's own description, falling back to the group's.
      header: 'Description',
      className: 'rg-table__cell--wrap',
      render: ({ asset }) => asset.description ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'size',
      header: 'Size',
      align: 'right',
      render: ({ asset }) => <FileSize bytes={asset.size} />,
    },
    {
      key: 'digest',
      header: 'Asset digest',
      render: ({ asset }) => <DigestChip digest={asset.digest} />,
    },
    {
      key: 'serving_modes',
      header: 'Serving modes',
      render: ({ asset }) => <ServingModeBadges modes={asset.serving_modes} />,
    },
    {
      key: 'downloads',
      header: 'Downloads',
      align: 'right',
      render: ({ archive }) =>
        archive ? archive.download_count : <span className="rg-muted">NA</span>,
    },
    {
      key: 'links',
      header: 'Links',
      render: ({ asset, archive }) => (
        <span className="flex flex-wrap items-center gap-2">
          {canDownloadArchives && archive && (
            <DownloadLink
              className="rg-btn rg-btn--sm"
              href={archiveDownloadUrl(client, archive.digest)}
            >
              Archive
            </DownloadLink>
          )}
          {asset.digest && (
            <ScopedLink className="rg-link text-sm" to={`/assets/${asset.digest}`}>
              Files
            </ScopedLink>
          )}
          {asset.recipe_id !== null && asset.recipe_id !== undefined && (
            <ScopedLink className="rg-link text-sm" to={`/recipes/${asset.recipe_id}`}>
              Recipe
            </ScopedLink>
          )}
          {asset.digest && (
            <DeleteAssetButton
              compact
              digest={asset.digest}
              registryPath={`${asset.asset_group_name ?? ''}:${asset.name}`}
              seekKeys={asset.seek_keys}
            />
          )}
          {/* Renders nothing unless the bridge is connected. */}
          {asset.genome_digest && asset.asset_group_name && (
            <BridgePullButton
              compact
              genomeDigest={asset.genome_digest}
              assetGroupName={asset.asset_group_name}
              assetName={asset.name}
            />
          )}
        </span>
      ),
    },
  ];

  return (
    <DataTable
      caption="Genome assets"
      columns={columns}
      rows={rows}
      rowKey={({ asset }, index) => asset.digest ?? `${asset.name}-${index}`}
      loading={loading}
      error={error}
      empty={empty}
      onRetry={onRetry}
    />
  );
}
