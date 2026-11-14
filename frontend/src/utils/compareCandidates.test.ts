import { describe, expect, it } from 'vitest';
import { rankCompareCandidates, sameSpecies } from './compareCandidates';
import type { CandidateGenome } from './compareCandidates';

const genome = (
  over: Partial<CandidateGenome> & { digest: string },
): CandidateGenome => ({
  aliases: [over.digest],
  species_name: 'Homo sapiens',
  taxon_id: 9606,
  assembly_source: null,
  assembly_accession: null,
  ...over,
});

const self = genome({
  digest: 'hg38',
  assembly_accession: 'GCA_000001405.15',
  assembly_source: 'ncbi',
});

describe('sameSpecies', () => {
  it('matches on taxon id, and falls back to the species name', () => {
    expect(sameSpecies(self, genome({ digest: 'x', species_name: 'other' }))).toBe(
      true,
    );
    expect(sameSpecies(self, genome({ digest: 'x', taxon_id: 10090 }))).toBe(false);
    expect(
      sameSpecies(
        genome({ digest: 'a', taxon_id: null }),
        genome({ digest: 'b', taxon_id: null, species_name: 'homo sapiens ' }),
      ),
    ).toBe(true);
    expect(
      sameSpecies(genome({ digest: 'a', taxon_id: null, species_name: null }), self),
    ).toBe(false);
  });
});

describe('rankCompareCandidates', () => {
  it('drops itself and other species, and puts the likeliest matches first', () => {
    const ranked = rankCompareCandidates(self, [
      self,
      genome({
        digest: 'mm10',
        taxon_id: 10090,
        assembly_accession: 'GCA_000001405.15',
      }),
      genome({ digest: 'alu', assembly_source: 'GenBank' }),
      genome({ digest: 'ucsc', assembly_source: 'UCSC' }),
      genome({
        digest: 'p14',
        assembly_accession: 'GCF_000001405.40',
        assembly_source: 'NCBI',
      }),
      genome({ digest: 'ensembl-like', assembly_source: 'NCBI' }),
      genome({ digest: 'full', assembly_accession: 'GCA_000001405.15' }),
    ]);
    expect(ranked.map((c) => c.digest)).toEqual([
      'full', // same accession
      'p14', // same assembly, other version
      'ensembl-like', // same source
      'alu',
      'ucsc',
    ]);
  });

  it('labels each choice with alias and accession or source', () => {
    const [first, second] = rankCompareCandidates(self, [
      genome({
        digest: 'a',
        aliases: ['hg38_full'],
        assembly_accession: 'GCA_000001405.15',
      }),
      genome({ digest: 'b', aliases: ['hg38_ucsc'], assembly_source: 'UCSC' }),
    ]);
    expect(first.label).toBe('hg38_full · GCA_000001405.15');
    expect(second.label).toBe('hg38_ucsc · UCSC');
  });
});
