import { Link } from 'react-router-dom';
import { DataTable } from '../common/DataTable';
import { ScopedLink } from '../common/ScopedLink';
import { speciesGenomesHref } from '../../utils/species';
import type { ReactNode } from 'react';
import type { Column } from '../common/DataTable';
import type { SpeciesRow } from '../../utils/species';

export interface SpeciesTableProps {
  species: SpeciesRow[] | undefined;
  loading?: boolean;
  error?: unknown;
  empty?: ReactNode;
  onRetry?: () => void;
  /** Render the "Show in tree" column. Off inside `/local`, which has no tree. */
  showTreeLink?: boolean;
}

export function SpeciesTable({
  species,
  loading,
  error,
  empty,
  onRetry,
  showTreeLink,
}: SpeciesTableProps) {
  const columns: Array<Column<SpeciesRow>> = [
    {
      key: 'species',
      header: 'Species',
      render: (row) => {
        const href = speciesGenomesHref(row);
        // Genomes with no species still count, but cannot be linked.
        if (!href) return <span className="rg-muted">Unspecified</span>;
        return (
          <ScopedLink className="rg-link font-medium" to={href}>
            {row.speciesName}
          </ScopedLink>
        );
      },
    },
    {
      key: 'common_name',
      header: 'Common name',
      render: (row) =>
        row.commonNames.length ? (
          row.commonNames.join(', ')
        ) : (
          <span className="rg-muted">NA</span>
        ),
    },
    {
      key: 'taxon_id',
      header: 'Taxon ID',
      // Plain text: the list response has no taxon URI to link to.
      render: (row) =>
        row.taxonIds.length ? row.taxonIds.join(', ') : <span className="rg-muted">NA</span>,
    },
    { key: 'genomes', header: 'Genomes', align: 'right', render: (row) => row.genomes },
    {
      key: 'assets',
      header: 'Assets',
      align: 'right',
      // Zero is styled as absence.
      render: (row) =>
        row.assets > 0 ? (
          <span className="font-semibold">{row.assets}</span>
        ) : (
          <span className="rg-muted">0</span>
        ),
    },
  ];

  if (showTreeLink) {
    columns.push({
      key: 'tree',
      header: 'Tree',
      // A plain Link, not a ScopedLink: there is no /local/tree.
      render: (row) =>
        row.key === '' ? (
          <span className="rg-muted">NA</span>
        ) : (
          <Link className="rg-link" to={`/tree?species=${encodeURIComponent(row.speciesName)}`}>
            Show in tree
          </Link>
        ),
    });
  }

  return (
    <DataTable
      caption="Species"
      columns={columns}
      rows={species}
      rowKey={(row) => row.key}
      loading={loading}
      error={error}
      empty={empty}
      onRetry={onRetry}
    />
  );
}
