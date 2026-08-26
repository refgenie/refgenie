import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { Badge } from '../common/Badge';
import { ActionBar } from '../common/ActionBar';
import type { Column } from '../common/DataTable';
import type { RemoteGenome } from '../../types/api';

export interface RemoteGenomeTableProps {
  genomes: RemoteGenome[] | undefined;
  /** Digests already present locally, for the "already local?" indicator. */
  localDigests: Set<string>;
  expanded: string | undefined;
  onToggle: (digest: string) => void;
  loading?: boolean;
  error?: unknown;
}

export function RemoteGenomeTable({
  genomes,
  localDigests,
  expanded,
  onToggle,
  loading,
  error,
}: RemoteGenomeTableProps) {
  const columns: Array<Column<RemoteGenome>> = [
    {
      key: 'aliases',
      header: 'Genome',
      render: (genome) => (
        <button
          type="button"
          className="rg-btn rg-btn--bare"
          onClick={() => onToggle(genome.genome_digest)}
          aria-expanded={expanded === genome.genome_digest}
        >
          {genome.aliases[0] ?? genome.genome_digest}
        </button>
      ),
    },
    {
      key: 'digest',
      header: 'Digest',
      render: (genome) => <DigestChip digest={genome.genome_digest} />,
    },
    {
      key: 'description',
      header: 'Description',
      render: (genome) => genome.description ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'local',
      header: 'Local',
      render: (genome) =>
        localDigests.has(genome.genome_digest) ? (
          <Badge variant="local">local</Badge>
        ) : (
          <span className="rg-muted">not local</span>
        ),
    },
    {
      key: 'actions',
      header: '',
      render: (genome) => (
        <ActionBar slot="remote-genome" context={{ digest: genome.genome_digest }} />
      ),
    },
  ];

  return (
    <DataTable
      caption="Remote genomes"
      columns={columns}
      rows={genomes}
      rowKey={(genome) => `${genome.server_url}-${genome.genome_digest}`}
      loading={loading}
      error={error}
    />
  );
}
