/** A VGP-style collection and three genomes, for the collection page tests. */

import type { GenomeDetailResponse, GenomeResponse } from '../types/api';
import type { CollectionConfig } from '../types/ui';

export const vgpCollection: CollectionConfig = {
  id: 'vgp',
  name: 'Vertebrate Genomes Project assemblies',
  short_name: 'VGP',
  blurb: 'Reference-quality vertebrate genome assemblies. Listed by the UCSC VGP hub.',
  note: 'An independent refgenie view of publicly available assemblies.',
  accent_color: '#2b7a78',
  logo_url: null,
  stores: ['vgp'],
  links: [
    {
      label: 'Vertebrate Genomes Project',
      url: 'https://vertebrategenomesproject.org/',
    },
    { label: 'GenomeArk', url: 'https://www.genomeark.org/' },
  ],
  genome_links: [
    {
      label: 'NCBI',
      template: 'https://www.ncbi.nlm.nih.gov/datasets/genome/{assembly_accession}/',
    },
    {
      label: 'GenomeArk',
      template: 'https://www.genomeark.org/genomeark-all/{species_underscored}.html',
    },
  ],
  columns: ['assembly_accession', 'haplotype', 'assembly_level'],
};

function genome(
  overrides: Partial<GenomeResponse> & { digest: string },
): GenomeResponse {
  return {
    aliases: [],
    description: null,
    asset_count: 0,
    species_name: null,
    common_name: null,
    synonyms: null,
    taxon_id: null,
    assembly_source: 'VGP',
    assembly_accession: null,
    assembly_level: null,
    species: null,
    ...overrides,
  };
}

export const vgpGenomes: GenomeResponse[] = [
  genome({
    digest: 'vgpdigest00000000000000000000001',
    aliases: ['fAciRut3.2_maternal_haplotype'],
    species_name: 'Acipenser ruthenus',
    assembly_accession: 'GCF_902713425.1',
    assembly_level: 'Chromosome',
    asset_count: 2,
  }),
  genome({
    digest: 'vgpdigest00000000000000000000002',
    aliases: ['aBomVar4.paternal.1'],
    species_name: 'Bombina variegata',
    assembly_accession: 'GCA_905165185.1',
    assembly_level: 'Scaffold',
  }),
  genome({
    digest: 'vgpdigest00000000000000000000003',
    aliases: ['fLatCha1.pri'],
    species_name: 'Acipenser ruthenus',
    assembly_accession: 'GCA_000225785.1',
    assembly_level: 'Scaffold',
  }),
];

export const vgpDetail: GenomeDetailResponse = {
  digest: vgpGenomes[0].digest,
  description: null,
  species_name: 'Acipenser ruthenus',
  common_name: null,
  synonyms: null,
  taxon_id: 7906,
  assembly_source: 'VGP',
  assembly_accession: 'GCF_902713425.1',
  assembly_level: 'Chromosome',
  remote_url: null,
  taxon_uri: null,
  species: null,
  fhr: { genome: 'fAciRut3.2_maternal_haplotype', version: '2' },
  refget_store_url: null,
};
