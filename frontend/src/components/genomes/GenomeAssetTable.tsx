import { Link } from 'react-router-dom';
import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { FileSize } from '../common/FileSize';
import { DeleteAssetButton } from '../actions/DeleteAssetButton';
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
 * The genome's assets, grouped by asset group.
 *
 * Driven from `/assets`, not `/archives`, and LEFT-joined to the archive
 * record: file-mode assets (e.g. `fasta_index` with `serving_modes: ["file"]`,
 * which have no tarball) must still get a row.
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
          <Link className="rg-link font-medium" to={`/assets/${asset.digest}`}>
            {asset.asset_group_name ? `${asset.asset_group_name}:` : ''}
            {asset.name}
          </Link>
        ) : (
          <span>{asset.name}</span>
        ),
    },
    {
      key: 'description',
      // The Jinja page labelled this "asset description" while showing the
      // GROUP's description. Show the asset's own, falling back to the group's.
      header: 'Description',
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
            <a className="rg-btn rg-btn--sm" href={archiveDownloadUrl(client, archive.digest)}>
              Archive
            </a>
          )}
          {asset.digest && (
            <Link className="rg-link text-sm" to={`/assets/${asset.digest}`}>
              Files
            </Link>
          )}
          {/* The Jinja page emitted /recipes/None when recipe_id was null. */}
          {asset.recipe_id !== null && asset.recipe_id !== undefined && (
            <Link className="rg-link text-sm" to={`/recipes/${asset.recipe_id}`}>
              Recipe
            </Link>
          )}
          {asset.digest && (
            <DeleteAssetButton
              compact
              digest={asset.digest}
              registryPath={`${asset.asset_group_name ?? ''}:${asset.name}`}
              seekKeys={asset.seek_keys}
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
