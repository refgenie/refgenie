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
  /** Built-in column keys to leave out, when extraColumns cover the same data. */
  hideColumns?: string[];
}

/** How many other aliases show under the genome name before "+N more". */
const OTHER_ALIASES_SHOWN = 2;

export function GenomeTable({
  genomes,
  loading,
  error,
  empty,
  onRetry,
  query,
  extraColumns,
  hideColumns,
}: GenomeTableProps) {
  const columns: Array<Column<GenomeResponse>> = [
    {
      key: 'aliases',
      header: 'Genome',
      className: 'rg-table__cell--wrap',
      // A real <Link> so keyboard and middle-click work. Other aliases sit on
      // a second line rather than in their own wide column.
      render: (genome) => {
        const shown = preferredAlias(genome.aliases, query);
        // Exclude whichever alias the primary cell chose, not blindly index 0.
        const rest = genome.aliases.filter((alias) => alias !== shown);
        const extra = rest.length - OTHER_ALIASES_SHOWN;
        return (
          <>
            <ScopedLink className="rg-link font-medium" to={`/genomes/${genome.digest}`}>
              {shown ?? formatDigest(genome.digest)}
            </ScopedLink>
            {rest.length > 0 && (
              <span className="rg-table__sub" title={rest.join('\n')}>
                {rest.slice(0, OTHER_ALIASES_SHOWN).join(', ')}
                {extra > 0 && `, +${extra} more`}
              </span>
            )}
          </>
        );
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
      className: 'rg-table__cell--wrap',
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
              {genome.species_name}
              <span className="rg-table__sub">{genome.common_name}</span>
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
      header: 'Built assets',
      align: 'right',
      // Zero is muted: most genomes have nothing built yet.
      render: (genome) =>
        genome.asset_count > 0 ? (
          <span className="font-semibold">{genome.asset_count}</span>
        ) : (
          <span className="rg-muted">0</span>
        ),
    },
  ];

  const shownColumns = hideColumns
    ? columns.filter((column) => !hideColumns.includes(column.key))
    : columns;

  return (
    <DataTable
      caption="Genomes"
      columns={extraColumns ? [...shownColumns, ...extraColumns] : shownColumns}
      rows={genomes}
      rowKey={(genome) => genome.digest}
      loading={loading}
      error={error}
      empty={empty}
      onRetry={onRetry}
    />
  );
}
