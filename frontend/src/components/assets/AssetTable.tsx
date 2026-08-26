import { Link } from 'react-router-dom';
import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { FileSize } from '../common/FileSize';
import { ServingModeBadges } from './ServingModeBadges';
import type { Column } from '../common/DataTable';
import type { AssetResponse } from '../../types/api';

export interface AssetTableProps {
  assets: AssetResponse[] | undefined;
  loading?: boolean;
  error?: unknown;
  empty?: React.ReactNode;
  onRetry?: () => void;
  /** Appended after the built-in columns (default-asset radio, delete, …). */
  extraColumns?: Array<Column<AssetResponse>>;
}

/** Flat asset list, used by the /assets index. */
export function AssetTable({
  assets,
  loading,
  error,
  empty,
  onRetry,
  extraColumns,
}: AssetTableProps) {
  const columns: Array<Column<AssetResponse>> = [
    {
      key: 'name',
      header: 'Asset',
      render: (asset) =>
        asset.digest ? (
          <Link className="rg-link" to={`/assets/${asset.digest}`}>
            {asset.asset_group_name ? `${asset.asset_group_name}:` : ''}
            {asset.name}
          </Link>
        ) : (
          <span>{asset.name}</span>
        ),
    },
    {
      key: 'genome',
      header: 'Genome',
      render: (asset) =>
        asset.genome_digest ? (
          <Link className="rg-link" to={`/genomes/${asset.genome_digest}`}>
            <DigestChip digest={asset.genome_digest} copyable={false} />
          </Link>
        ) : (
          <span className="rg-muted">NA</span>
        ),
    },
    {
      key: 'asset_class',
      header: 'Class',
      render: (asset) => asset.asset_class_name ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'description',
      header: 'Description',
      render: (asset) => asset.description ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'serving_modes',
      header: 'Serving modes',
      render: (asset) => <ServingModeBadges modes={asset.serving_modes} />,
    },
    {
      key: 'size',
      header: 'Size',
      align: 'right',
      render: (asset) => <FileSize bytes={asset.size} />,
    },
    {
      key: 'digest',
      header: 'Digest',
      render: (asset) => <DigestChip digest={asset.digest} />,
    },
  ];

  return (
    <DataTable
      caption="Assets"
      columns={extraColumns ? [...columns, ...extraColumns] : columns}
      rows={assets}
      rowKey={(asset, index) => asset.digest ?? `${asset.name}-${index}`}
      loading={loading}
      error={error}
      empty={empty}
      onRetry={onRetry}
    />
  );
}
