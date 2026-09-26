import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { FileSize } from '../common/FileSize';
import { Badge } from '../common/Badge';
import type { Column } from '../common/DataTable';
import type { StagedAssetPublic } from '../../types/api';

export interface StagedAssetTableProps {
  staged: StagedAssetPublic[] | undefined;
  loading?: boolean;
  error?: unknown;
}

export function StagedAssetTable({ staged, loading, error }: StagedAssetTableProps) {
  const columns: Array<Column<StagedAssetPublic>> = [
    {
      key: 'mode',
      header: 'Mode',
      render: (record) => <Badge variant={record.mode}>{record.mode}</Badge>,
    },
    {
      key: 'asset_digest',
      header: 'Asset digest',
      render: (record) => <DigestChip digest={record.asset_digest} />,
    },
    {
      key: 'tarball_digest',
      header: 'Tarball digest',
      render: (record) => <DigestChip digest={record.tarball_digest} copyable={false} />,
    },
    {
      key: 'tarball_size',
      header: 'Tarball size',
      align: 'right',
      render: (record) => <FileSize bytes={record.tarball_size} />,
    },
    {
      key: 'download_count',
      header: 'Downloads',
      align: 'right',
      render: (record) => record.download_count,
    },
    {
      key: 'contents',
      header: 'Contents',
      render: (record) =>
        record.directory_contents?.length ? (
          <span className="text-xs rg-muted">{record.directory_contents.join(', ')}</span>
        ) : (
          <span className="rg-muted">NA</span>
        ),
    },
  ];

  return (
    <DataTable
      caption="Staged assets"
      columns={columns}
      rows={staged}
      rowKey={(record) => `${record.asset_digest}-${record.mode}`}
      loading={loading}
      error={error}
    />
  );
}
