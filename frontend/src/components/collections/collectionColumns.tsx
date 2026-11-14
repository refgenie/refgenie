import type { Column } from '../common/DataTable';
import { ExternalLink } from '../common/ExternalLink';
import { expandGenomeLink, usesPlaceholder } from '../../utils/collections';
import { haplotypeOf } from '../../utils/haplotype';
import { preferredAlias } from '../../utils/genomes';
import { formatDigest } from '../../utils/format';
import type { GenomeResponse } from '../../types/api';
import type { CollectionColumn, CollectionConfig } from '../../types/ui';

const NA = <span className="rg-muted">NA</span>;

function accessionCell(genome: GenomeResponse, collection: CollectionConfig) {
  if (!genome.assembly_accession) return null;
  const template = collection.genome_links.find((l) =>
    usesPlaceholder(l.template, 'assembly_accession'),
  );
  const href = template && expandGenomeLink(template.template, genome);
  return href ? (
    <ExternalLink href={href} hint={`${template.label}, opens in a new tab`}>
      {genome.assembly_accession}
    </ExternalLink>
  ) : (
    genome.assembly_accession
  );
}

/**
 * Accession and source share one "Assembly" column, standing in for the
 * built-in one: the accession (linked when the collection has a template for
 * it) with the source beneath.
 */
function assemblyColumn(collection: CollectionConfig): Column<GenomeResponse> {
  return {
    key: 'collection_assembly',
    header: 'Assembly',
    render: (genome) => {
      const accession = accessionCell(genome, collection);
      const source = genome.assembly_source;
      if (!accession) return source ?? NA;
      return (
        <>
          {accession}
          {source && <span className="rg-table__sub">{source}</span>}
        </>
      );
    },
  };
}

function column(
  key: Exclude<CollectionColumn, 'assembly_accession' | 'assembly_source'>,
): Column<GenomeResponse> {
  switch (key) {
    case 'haplotype':
      return {
        key,
        header: 'Haplotype',
        render: (genome) => {
          const haplotype = haplotypeOf(genome);
          return haplotype ? (
            <span className="rg-badge rg-badge--default rg-collection-haplotype">
              {haplotype}
            </span>
          ) : (
            NA
          );
        },
      };
    case 'assembly_level':
      return {
        key,
        header: 'Assembly level',
        render: (genome) => genome.assembly_level ?? NA,
      };
  }
}

/** Built-in GenomeTable columns a collection table leaves out: its own columns cover them. */
export const COLLECTION_HIDDEN_COLUMNS = ['assembly', 'description'];

/** The collection's configured columns, then a Details button, for `GenomeTable.extraColumns`. */
export function collectionColumns(
  collection: CollectionConfig,
  onDetails: (genome: GenomeResponse) => void,
): Array<Column<GenomeResponse>> {
  // The accession and source keys both fold into the one Assembly column.
  const rest = collection.columns.filter(
    (key): key is Exclude<CollectionColumn, 'assembly_accession' | 'assembly_source'> =>
      key !== 'assembly_accession' && key !== 'assembly_source',
  );
  return [
    assemblyColumn(collection),
    ...rest.map(column),
    {
      key: 'details',
      header: 'Details',
      render: (genome) => (
        <button
          type="button"
          className="rg-btn rg-btn--sm"
          aria-label={`Details for ${preferredAlias(genome.aliases) ?? formatDigest(genome.digest)}`}
          onClick={() => onDetails(genome)}
        >
          Details
        </button>
      ),
    },
  ];
}
