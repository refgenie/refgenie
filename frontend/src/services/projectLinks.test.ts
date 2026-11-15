import { describe, expect, it } from 'vitest';
import { genomeGuide } from './projectLinks';

describe('genomeGuide', () => {
  it('finds the guide for an alias, a species or a search', () => {
    expect(genomeGuide(['hs38DH', 'GRCh38-full-decoy-hla'])?.href).toMatch(
      /choosing_hg38\/$/,
    );
    expect(genomeGuide(['mm39_ucsc'])?.title).toBe('Which mm39 should I use?');
    expect(genomeGuide(['ASM1', 'Saccharomyces cerevisiae'])?.href).toMatch(
      /choosing_yeast\/$/,
    );
    expect(genomeGuide([' Yeast '])?.title).toBe('Which yeast genome should I use?');
  });

  it('has none for other genomes', () => {
    expect(genomeGuide(['hg19', 'Homo sapiens', null, undefined, ''])).toBeUndefined();
    expect(genomeGuide(['mm10'])).toBeUndefined();
  });
});
