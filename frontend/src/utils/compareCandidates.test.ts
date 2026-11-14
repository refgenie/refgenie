import { describe, expect, it } from 'vitest';
import { nameStem, rankCompareCandidates, sameSpecies } from './compareCandidates';
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

  it('puts genomes that share a name stem first, by name', () => {
    const hg38 = genome({
      digest: 'hg38',
      aliases: ['hg38', 'GRCh38-fasta-no-alt-analysis'],
      assembly_accession: 'GCA_000001405.15',
      assembly_source: 'NCBI',
    });
    const ranked = rankCompareCandidates(hg38, [
      genome({ digest: 'a', aliases: ['hg19'], assembly_accession: 'GCF_000001405.25' }),
      genome({ digest: 'b', aliases: ['hg38_primary'], assembly_source: 'UCSC' }),
      genome({
        digest: 'c',
        aliases: ['hg38_full_analysis'],
        assembly_accession: 'GCA_000001405.15',
      }),
      genome({ digest: 'd', aliases: ['GRCh38-81-ensembl', 'homo_sapiens_ENSEMBL_release_81'] }),
      genome({ digest: 'e', aliases: ['hg38-ddbj', 'hg38_initial'] }),
      // Its name is hg19's; a stray GRCh38 alias does not pull it up.
      genome({ digest: 'f', aliases: ['hg19_b37', 'GRCh38-b37-mislabeled'] }),
    ]);
    expect(ranked.map((c) => c.name)).toEqual([
      'GRCh38-81-ensembl', // shares the GRCh38 stem with an alias of hg38
      'hg38_full_analysis',
      'hg38_primary',
      'hg38-ddbj',
      'hg19', // not the family: same assembly accession, other version
      'hg19_b37',
    ]);
    // Every alias rides along, for type-to-filter.
    expect(ranked[0].aliases).toEqual(['GRCh38-81-ensembl', 'homo_sapiens_ENSEMBL_release_81']);
  });
});

describe('nameStem', () => {
  it('is the leading letters and digits, lowercased', () => {
    expect(nameStem('hg38_primary')).toBe('hg38');
    expect(nameStem('GRCh38.p14-fasta')).toBe('grch38');
    expect(nameStem('rCRSd')).toBeNull();
  });
});
