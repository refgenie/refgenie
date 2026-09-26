import { Link } from 'react-router-dom';
import { DataTable } from '../common/DataTable';
import { ServingModeBadges } from '../assets/ServingModeBadges';
import type { Column } from '../common/DataTable';
import type { AssetClassPublic } from '../../types/api';

export interface AssetClassTableProps {
  assetClasses: AssetClassPublic[] | undefined;
  loading?: boolean;
  error?: unknown;
  empty?: React.ReactNode;
  onRetry?: () => void;
}

export function AssetClassTable({
  assetClasses,
  loading,
  error,
  empty,
  onRetry,
}: AssetClassTableProps) {
  const columns: Array<Column<AssetClassPublic>> = [
    {
      key: 'name',
      header: 'Name',
      render: (item) =>
        item.id !== null ? (
          <Link className="rg-link font-medium" to={`/asset-classes/${item.id}`}>
            {item.name}
          </Link>
        ) : (
          <span>{item.name}</span>
        ),
    },
    { key: 'version', header: 'Version', render: (item) => item.version },
    {
      key: 'description',
      header: 'Description',
      render: (item) => item.description ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'serving_modes',
      header: 'Serving modes',
      render: (item) => <ServingModeBadges modes={item.serving_modes} />,
    },
  ];

  return (
    <DataTable
      caption="Asset classes"
      columns={columns}
      rows={assetClasses}
      rowKey={(item, index) => String(item.id ?? `${item.name}-${index}`)}
      loading={loading}
      error={error}
      empty={empty}
      onRetry={onRetry}
    />
  );
}
