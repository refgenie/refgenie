import { BaseModal } from '../common/BaseModal';
import { DescriptionList } from '../common/DescriptionList';
import { ExternalLink } from '../common/ExternalLink';
import { ScopedLink } from '../common/ScopedLink';
import { LoadingState, ErrorState } from '../common/states';
import { FhrPanel } from '../genomes/FhrPanel';
import { useGenome } from '../../hooks/queries/useGenomes';
import { expandGenomeLink } from '../../utils/collections';
import { preferredAlias } from '../../utils/genomes';
import { formatDigest } from '../../utils/format';
import { haplotypeOf } from '../../utils/haplotype';
import type { GenomeResponse } from '../../types/api';
import type { CollectionConfig } from '../../types/ui';

export interface CollectionGenomeDetailProps {
  collection: CollectionConfig;
  /** The listing row; metadata is fetched from the detail endpoint. */
  genome: GenomeResponse;
  onClose: () => void;
}

/** One genome of a collection: identifiers, external links and metadata. */
export function CollectionGenomeDetail({
  collection,
  genome,
  onClose,
}: CollectionGenomeDetailProps) {
  const detail = useGenome(genome.digest);
  const title = preferredAlias(genome.aliases) ?? formatDigest(genome.digest);
  const links = collection.genome_links
    .map((link) => ({
      label: link.label,
      href: expandGenomeLink(link.template, genome),
    }))
    .filter((link): link is { label: string; href: string } => link.href !== null);
  const fhr = detail.data?.fhr;
  const hasVitalStats = !!fhr?.vitalStats && Object.keys(fhr.vitalStats).length > 0;

  return (
    <BaseModal isOpen onClose={onClose} title={title} size="lg">
      <div className="flex flex-col gap-4">
        <DescriptionList
          items={[
            { term: 'Description', value: genome.description ?? '' },
            { term: 'Species', value: genome.species_name ?? '' },
            { term: 'Accession', value: genome.assembly_accession ?? '' },
            { term: 'Haplotype', value: haplotypeOf(genome) ?? '' },
            { term: 'Assembly level', value: genome.assembly_level ?? '' },
            { term: 'Source', value: genome.assembly_source ?? '' },
            {
              term: 'Links',
              value: links.length ? (
                <span className="flex flex-wrap gap-3">
                  {links.map((link) => (
                    <ExternalLink key={link.label} href={link.href}>
                      {link.label}
                    </ExternalLink>
                  ))}
                </span>
              ) : (
                ''
              ),
            },
          ]}
        />
        <ScopedLink className="rg-link" to={`/genomes/${genome.digest}`}>
          Open the full genome page
        </ScopedLink>
        {detail.isPending ? (
          <LoadingState label="Loading metadata" />
        ) : detail.error ? (
          <ErrorState error={detail.error} />
        ) : (
          <>
            {!hasVitalStats && (
              <p className="rg-muted text-sm">
                No quality metrics have been published for this assembly.
              </p>
            )}
            {fhr && <FhrPanel fhr={fhr} />}
          </>
        )}
      </div>
    </BaseModal>
  );
}
