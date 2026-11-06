import { DataTable } from '../common/DataTable';
import { ScopedLink } from '../common/ScopedLink';
import { DigestChip } from '../common/DigestChip';
import { preferredAlias } from '../../utils/genomes';
import { formatDigest } from '../../utils/format';
import type { Column } from '../common/DataTable';
import type { GenomeResponse } from '../../types/api';

export interface GenomeTableProps {
  genomes: GenomeResponse[] | undefined;
  loading?: boolean;
  error?: unknown;
  empty?: React.ReactNode;
  onRetry?: () => void;
  /** The live search term; each row is labelled with the alias that matched. */
  query?: string;
  /** Appended after the built-in columns (the bridge's local-presence badge). */
  extraColumns?: Array<Column<GenomeResponse>>;
}

export function GenomeTable({
  genomes,
  loading,
  error,
  empty,
  onRetry,
  query,
  extraColumns,
}: GenomeTableProps) {
  const columns: Array<Column<GenomeResponse>> = [
    {
      key: 'aliases',
      header: 'Genome',
      // A real <Link> so keyboard and middle-click work.
      render: (genome) => (
        <ScopedLink className="rg-link font-medium" to={`/genomes/${genome.digest}`}>
          {preferredAlias(genome.aliases, query) ?? formatDigest(genome.digest)}
        </ScopedLink>
      ),
    },
    {
      key: 'other_aliases',
      header: 'Other aliases',
      // Exclude whichever alias the primary cell chose, not blindly index 0.
      render: (genome) => {
        const shown = preferredAlias(genome.aliases, query);
        const rest = genome.aliases.filter((alias) => alias !== shown);
        return rest.length ? rest.join(', ') : <span className="rg-muted">NA</span>;
      },
    },
    {
      key: 'digest',
      header: 'Digest',
      render: (genome) => <DigestChip digest={genome.digest} />,
    },
    {
      key: 'description',
      header: 'Description',
      render: (genome) => genome.description ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'species',
      header: 'Species',
      // One common name beside the scientific name.
      render: (genome) => {
        if (genome.species_name && genome.common_name) {
          return (
            <>
              {genome.species_name}{' '}
              <span className="rg-muted">{genome.common_name}</span>
            </>
          );
        }
        return (
          genome.species_name ??
          genome.common_name ?? <span className="rg-muted">NA</span>
        );
      },
    },
    {
      key: 'assembly',
      header: 'Assembly',
      render: (genome) => {
        const parts = [genome.assembly_source, genome.assembly_accession].filter(
          Boolean,
        );
        return parts.length ? parts.join(' · ') : <span className="rg-muted">NA</span>;
      },
    },
    {
      key: 'assets',
      header: 'Assets',
      align: 'right',
      // Most genomes have nothing built, so zero is styled as absence.
      render: (genome) =>
        genome.asset_count > 0 ? (
          <span className="font-semibold">{genome.asset_count}</span>
        ) : (
          <span className="rg-muted">0</span>
        ),
    },
  ];

  return (
    <DataTable
      caption="Genomes"
      columns={extraColumns ? [...columns, ...extraColumns] : columns}
      rows={genomes}
      rowKey={(genome) => genome.digest}
      loading={loading}
      error={error}
      empty={empty}
      onRetry={onRetry}
    />
  );
}
