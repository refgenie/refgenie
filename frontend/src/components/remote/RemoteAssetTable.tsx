import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { FileSize } from '../common/FileSize';
import { Badge } from '../common/Badge';
import { PullButton } from '../actions/PullButton';
import type { Column } from '../common/DataTable';
import type { RemoteAsset } from '../../types/api';

export interface RemoteAssetTableProps {
  assets: RemoteAsset[] | undefined;
  /** Asset digests already present locally. */
  localDigests: Set<string>;
  loading?: boolean;
  error?: unknown;
  /** Offered as "try another server" when more than one server is subscribed. */
  onPickServer?: () => void;
}

export function RemoteAssetTable({
  assets,
  localDigests,
  loading,
  error,
  onPickServer,
}: RemoteAssetTableProps) {
  const columns: Array<Column<RemoteAsset>> = [
    {
      key: 'name',
      header: 'Asset',
      render: (asset) => `${asset.asset_group_name}:${asset.asset_name}`,
    },
    {
      key: 'size',
      header: 'Archive size',
      align: 'right',
      render: (asset) => <FileSize bytes={asset.archive_size} />,
    },
    {
      key: 'digest',
      header: 'Asset digest',
      render: (asset) => <DigestChip digest={asset.asset_digest} />,
    },
    {
      key: 'local',
      header: 'Local',
      render: (asset) =>
        localDigests.has(asset.asset_digest) ? (
          <Badge variant="local">local</Badge>
        ) : (
          <span className="rg-muted">not local</span>
        ),
    },
    {
      key: 'actions',
      header: '',
      render: (asset) => (
        <PullButton
          serverUrl={asset.server_url}
          genomeDigest={asset.genome_digest}
          assetGroupName={asset.asset_group_name}
          assetName={asset.asset_name}
          assetDigest={asset.asset_digest}
          archiveDigest={asset.archive_digest}
          archiveSize={asset.archive_size}
          existsLocally={localDigests.has(asset.asset_digest)}
          onPickServer={onPickServer}
        />
      ),
    },
  ];

  return (
    <DataTable
      caption="Remote assets"
      columns={columns}
      rows={assets}
      rowKey={(asset) => `${asset.server_url}-${asset.asset_digest}`}
      loading={loading}
      error={error}
    />
  );
}
