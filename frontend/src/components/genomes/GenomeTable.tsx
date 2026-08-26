import { Link } from 'react-router-dom';
import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import type { Column } from '../common/DataTable';
import type { GenomeResponse } from '../../types/api';

export interface GenomeTableProps {
  genomes: GenomeResponse[] | undefined;
  loading?: boolean;
  error?: unknown;
  empty?: React.ReactNode;
  onRetry?: () => void;
}

function primaryLabel(genome: GenomeResponse): string {
  return genome.aliases[0] ?? genome.digest;
}

export function GenomeTable({ genomes, loading, error, empty, onRetry }: GenomeTableProps) {
  const columns: Array<Column<GenomeResponse>> = [
    {
      key: 'aliases',
      header: 'Genome',
      // A real <Link> in the primary cell so keyboard users and middle-click
      // work regardless of any whole-row click handler.
      render: (genome) => (
        <Link className="rg-link font-medium" to={`/genomes/${genome.digest}`}>
          {primaryLabel(genome)}
        </Link>
      ),
    },
    {
      key: 'other_aliases',
      header: 'Other aliases',
      render: (genome) =>
        genome.aliases.length > 1 ? (
          genome.aliases.slice(1).join(', ')
        ) : (
          <span className="rg-muted">NA</span>
        ),
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
      render: (genome) =>
        genome.species_name ?? genome.common_name ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'assembly',
      header: 'Assembly',
      render: (genome) => {
        const parts = [genome.assembly_source, genome.assembly_accession].filter(Boolean);
        return parts.length ? parts.join(' · ') : <span className="rg-muted">NA</span>;
      },
    },
    {
      key: 'assets',
      header: 'Assets',
      align: 'right',
      render: (genome) => genome.asset_count,
    },
    {
      key: 'fasta',
      header: 'fasta',
      // Shortcut to the default fasta asset, matching the Jinja index column.
      render: (genome) => (
        <Link className="rg-link" to={`/assets?genome_digest=${genome.digest}&q=fasta`}>
          fasta
        </Link>
      ),
    },
  ];

  return (
    <DataTable
      caption="Genomes"
      columns={columns}
      rows={genomes}
      rowKey={(genome) => genome.digest}
      loading={loading}
      error={error}
      empty={empty}
      onRetry={onRetry}
    />
  );
}
