import type { Column } from '../common/DataTable';
import { ExternalLink } from '../common/ExternalLink';
import { expandGenomeLink, usesPlaceholder } from '../../utils/collections';
import { haplotypeOf } from '../../utils/haplotype';
import { preferredAlias } from '../../utils/genomes';
import { formatDigest } from '../../utils/format';
import type { GenomeResponse } from '../../types/api';
import type { CollectionColumn, CollectionConfig } from '../../types/ui';

const NA = <span className="rg-muted">NA</span>;

function column(
  key: CollectionColumn,
  collection: CollectionConfig,
): Column<GenomeResponse> {
  switch (key) {
    case 'assembly_accession': {
      const template = collection.genome_links.find((l) =>
        usesPlaceholder(l.template, 'assembly_accession'),
      );
      return {
        key,
        header: 'Accession',
        render: (genome) => {
          if (!genome.assembly_accession) return NA;
          const href = template && expandGenomeLink(template.template, genome);
          return href ? (
            <ExternalLink href={href} hint={`${template.label}, opens in a new tab`}>
              {genome.assembly_accession}
            </ExternalLink>
          ) : (
            genome.assembly_accession
          );
        },
      };
    }
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
    case 'assembly_source':
      return {
        key,
        header: 'Source',
        render: (genome) => genome.assembly_source ?? NA,
      };
  }
}

/** The collection's configured columns, then a Details button, for `GenomeTable.extraColumns`. */
export function collectionColumns(
  collection: CollectionConfig,
  onDetails: (genome: GenomeResponse) => void,
): Array<Column<GenomeResponse>> {
  return [
    ...collection.columns.map((key) => column(key, collection)),
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
