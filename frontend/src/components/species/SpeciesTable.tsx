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
  /**
   * Render the "Show in tree" column. The tree lives on the page's own backend
   * and is capability-gated, so the caller decides; a tree link fired from
   * inside `/local` would silently cross backends.
   */
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
        // A genome with no species recorded still counts; it just cannot be
        // linked, because no query parameter expresses "species_name IS NULL".
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
      // Plain text, not a link: the server hands out a `taxon_uri` on the
      // genome DETAIL response, and inventing an NCBI URL here would be the UI
      // guessing at a resolver the backend never named.
      render: (row) =>
        row.taxonIds.length ? row.taxonIds.join(', ') : <span className="rg-muted">NA</span>,
    },
    { key: 'genomes', header: 'Genomes', align: 'right', render: (row) => row.genomes },
    {
      key: 'assets',
      header: 'Assets',
      align: 'right',
      // A zero is the answer to "can I download anything for this organism?",
      // so it reads as absence rather than as one more number in the column.
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
      // A plain Link, not a ScopedLink: there is no /local/tree, so this one
      // intentionally exits the browse scope — the same rule asset-class and
      // recipe links already follow.
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
