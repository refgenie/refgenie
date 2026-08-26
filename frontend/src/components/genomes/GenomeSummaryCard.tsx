import { DescriptionList } from '../common/DescriptionList';
import { DigestChip } from '../common/DigestChip';
import { ExternalLink } from '../common/ExternalLink';
import { accessionUrl, taxonUrl } from '../../utils/links';
import type { GenomeDetailResponse } from '../../types/api';

export interface GenomeSummaryCardProps {
  genome: GenomeDetailResponse;
  aliases: string[];
}

export function GenomeSummaryCard({ genome, aliases }: GenomeSummaryCardProps) {
  const taxon = taxonUrl(genome.taxon_uri, genome.taxon_id);
  const accession = accessionUrl(
    genome.assembly_accession,
    typeof genome.fhr?.accessionID?.url === 'string' ? genome.fhr.accessionID.url : null,
  );

  return (
    <DescriptionList
      items={[
        {
          term: 'Aliases',
          value: aliases.length ? <strong>{aliases.join(', ')}</strong> : '',
        },
        { term: 'Digest', value: <DigestChip digest={genome.digest} length={32} /> },
        { term: 'Description', value: genome.description ?? '' },
        { term: 'Scientific name', value: genome.species_name ?? '' },
        { term: 'Common name', value: genome.common_name ?? '' },
        {
          term: 'NCBI taxon',
          value:
            genome.taxon_id === null || genome.taxon_id === undefined ? (
              ''
            ) : taxon ? (
              <ExternalLink href={taxon}>{genome.taxon_id}</ExternalLink>
            ) : (
              String(genome.taxon_id)
            ),
        },
        { term: 'Assembly source', value: genome.assembly_source ?? '' },
        {
          term: 'Assembly accession',
          value: !genome.assembly_accession ? (
            ''
          ) : accession ? (
            <ExternalLink href={accession}>{genome.assembly_accession}</ExternalLink>
          ) : (
            genome.assembly_accession
          ),
        },
        { term: 'Assembly level', value: genome.assembly_level ?? '' },
        {
          term: 'Remote URL',
          value: genome.remote_url ? (
            <ExternalLink href={genome.remote_url}>{genome.remote_url}</ExternalLink>
          ) : (
            ''
          ),
        },
      ]}
    />
  );
}
