import { ScopedLink } from '../common/ScopedLink';
import { useAssetGroups } from '../../hooks/queries/useAssetGroups';
import { useGenomeIndex } from '../../hooks/queries/useGenomes';
import { useResetPastEndPage } from '../../hooks/useResetPastEndPage';
import { genomeLabel } from '../../utils/genomes';
import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { Pagination } from '../common/Pagination';
import { EmptyState } from '../common/states';
import type { Column } from '../common/DataTable';
import type { AssetGroupPublic } from '../../types/api';

export interface AssetClassGenomeTableProps {
  /** The route's numeric id, used to drop other versions of the same name. */
  assetClassId: number;
  /** The only thing `/asset_groups` can filter on. */
  assetClassName: string;
  /** 0-based. */
  page: number;
  pageSize: number;
  onPageChange: (page: number) => void;
}

/**
 * Every genome that has this asset class built, from
 * `/asset_groups?asset_class=`. Aliases come from the cached genome index.
 */
export function AssetClassGenomeTable({
  assetClassId,
  assetClassName,
  page,
  pageSize,
  onPageChange,
}: AssetClassGenomeTableProps) {
  const groups = useAssetGroups({ asset_class: assetClassName, page, pageSize });
  useResetPastEndPage(groups.error, page, onPageChange);
  const genomes = useGenomeIndex();
  const index = genomes.data;

  // The endpoint filters by class name only, so drop rows from other versions
  // (the total still counts them). It has no `sort`, so this sorts one page.
  const rows = (groups.data?.results ?? [])
    .filter((group) => group.asset_class_id === assetClassId)
    .sort((a, b) =>
      genomeLabel(index, a.genome_digest).localeCompare(genomeLabel(index, b.genome_digest)),
    );

  const columns: Array<Column<AssetGroupPublic>> = [
    {
      key: 'genome',
      header: 'Genome',
      render: (group) => (
        <ScopedLink className="rg-link font-medium" to={`/genomes/${group.genome_digest}`}>
          {genomeLabel(index, group.genome_digest)}
        </ScopedLink>
      ),
    },
    {
      key: 'species',
      header: 'Species',
      render: (group) => {
        const genome = index?.get(group.genome_digest);
        const label = genome?.species_name ?? genome?.common_name;
        return label ?? <span className="rg-muted">NA</span>;
      },
    },
    {
      key: 'assets',
      header: 'Assets',
      render: (group) =>
        group.id !== null ? (
          <ScopedLink className="rg-link" to={`/asset-groups/${group.id}`}>
            {group.name}
          </ScopedLink>
        ) : (
          <span>{group.name}</span>
        ),
    },
    {
      key: 'digest',
      header: 'Genome digest',
      render: (group) => <DigestChip digest={group.genome_digest} />,
    },
  ];

  return (
    <section>
      <h2 className="text-xl font-semibold mb-4">
        Genomes
        {groups.data && (
          <span className="rg-muted text-sm font-normal">
            {' '}
            ({groups.data.pagination.total} with{' '}
            <code className="rg-code rg-code--inline">{assetClassName}</code> built)
          </span>
        )}
      </h2>
      <DataTable
        caption="Genomes with this asset class"
        columns={columns}
        rows={rows}
        rowKey={(group, i) => String(group.id ?? `${group.genome_digest}-${i}`)}
        loading={groups.isPending}
        error={groups.error}
        onRetry={() => groups.refetch()}
        empty={<EmptyState message="No genome on this server has this asset class built." />}
      />
      <Pagination pagination={groups.data?.pagination} onPageChange={onPageChange} />
    </section>
  );
}
