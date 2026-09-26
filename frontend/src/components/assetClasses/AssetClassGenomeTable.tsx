import { Link } from 'react-router-dom';
import { useAssetGroups } from '../../hooks/queries/useAssetGroups';
import { useGenomeIndex } from '../../hooks/queries/useGenomes';
import { genomeLabel } from '../../utils/genomes';
import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { Pagination } from '../common/Pagination';
import { EmptyState } from '../common/states';
import type { Column } from '../common/DataTable';
import type { AssetGroupPublic } from '../../types/api';

export interface AssetClassGenomeTableProps {
  /** The route's numeric id — the exact class, used to guard the name filter. */
  assetClassId: number;
  /** The only thing `/asset_groups` can filter on. */
  assetClassName: string;
  offset: number;
  limit: number;
  onOffsetChange: (offset: number) => void;
}

/**
 * Every genome that has this asset class built.
 *
 * Driven from `/asset_groups?asset_class=`, because an AssetGroup IS the
 * (genome x asset class) pairing. The rows carry only `genome_digest`, so the
 * alias comes from the cached genome index; the row's own id is the deep link
 * to `/asset-groups/:id`, which lists that pair's assets with sizes.
 */
export function AssetClassGenomeTable({
  assetClassId,
  assetClassName,
  offset,
  limit,
  onOffsetChange,
}: AssetClassGenomeTableProps) {
  const groups = useAssetGroups({ asset_class: assetClassName, offset, limit });
  const genomes = useGenomeIndex();
  const index = genomes.data;

  // `/asset_groups?asset_class=` filters on AssetClass.NAME, and AssetClass is
  // unique on (name, version). Drop rows belonging to a sibling version.
  //
  // Consequence: `pagination.total` still counts every version, so a page can
  // render short. No asset class name is duplicated on the public server today,
  // so this is a correctness guard rather than a live condition.
  //
  // The sort is within the page only — `/asset_groups` has no `sort` parameter
  // and returns insertion order. Every class on the public server fits in one
  // 50-row page (the largest, `fasta`, has 26), so in practice this is a full
  // ordering; on a server where a class spills to a second page it is not.
  const rows = (groups.data?.items ?? [])
    .filter((group) => group.asset_class_id === assetClassId)
    .sort((a, b) =>
      genomeLabel(index, a.genome_digest).localeCompare(genomeLabel(index, b.genome_digest)),
    );

  const columns: Array<Column<AssetGroupPublic>> = [
    {
      key: 'genome',
      header: 'Genome',
      render: (group) => (
        <Link className="rg-link font-medium" to={`/genomes/${group.genome_digest}`}>
          {genomeLabel(index, group.genome_digest)}
        </Link>
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
      // The asset group page is the exact (genome x class) pair: it lists the
      // built assets with sizes and per-asset links, so each row deep-links
      // there.
      render: (group) =>
        group.id !== null ? (
          <Link className="rg-link" to={`/asset-groups/${group.id}`}>
            {group.name}
          </Link>
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
            — {groups.data.pagination.total} with{' '}
            <code className="rg-code rg-code--inline">{assetClassName}</code> built
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
      <Pagination pagination={groups.data?.pagination} onOffsetChange={onOffsetChange} />
    </section>
  );
}
