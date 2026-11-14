import { useAssetClasses } from '../hooks/queries/useAssetClasses';
import { ActionBar } from '../components/common/ActionBar';
import { ScopedLink } from '../components/common/ScopedLink';
import { SearchListPage } from '../components/common/SearchListPage';
import { ServingModeBadges } from '../components/assets/ServingModeBadges';
import { ASSET_CLASS_SEARCH_FIELDS } from '../services/resources/assetClasses';
import type { Column } from '../components/common/DataTable';
import type { AssetClassPublic } from '../types/api';

const COLUMNS: Array<Column<AssetClassPublic>> = [
  {
    key: 'name',
    header: 'Name',
    render: (item) =>
      item.id !== null ? (
        <ScopedLink className="rg-link font-medium" to={`/asset-classes/${item.id}`}>
          {item.name}
        </ScopedLink>
      ) : (
        <span>{item.name}</span>
      ),
  },
  { key: 'version', header: 'Version', render: (item) => item.version },
  {
    key: 'description',
    header: 'Description',
    className: 'rg-table__cell--wrap',
    render: (item) => item.description ?? <span className="rg-muted">NA</span>,
  },
  {
    key: 'serving_modes',
    header: 'Serving modes',
    render: (item) => <ServingModeBadges modes={item.serving_modes} />,
  },
];

export function AssetClassesPage() {
  return (
    <SearchListPage
      title="Asset classes"
      actions={<ActionBar slot="asset-class" />}
      lede={
        <>
          An asset class is the definition of one kind of genome file bundle, such as a
          FASTA, a bowtie2 index, or a set of annotations. It names the files that bundle must contain, so
          anything built to that class can be used the same way whichever genome it came from.
        </>
      }
      useList={useAssetClasses}
      searchFields={ASSET_CLASS_SEARCH_FIELDS}
      columns={COLUMNS}
      rowKey={(item, index) => String(item.id ?? `${item.name}-${index}`)}
      emptyMessage="No asset classes registered."
    />
  );
}
