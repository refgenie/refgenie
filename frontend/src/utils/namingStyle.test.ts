import { describe, expect, it } from 'vitest';
import { namingLabel, namingStyle } from './namingStyle';

const seqs = (names: string[]) =>
  names.map((name, i) => ({ name, length: 1_000_000 - i * 1000 }));

describe('namingStyle', () => {
  it('reads UCSC names, with the mitochondrion as an example', () => {
    const naming = namingStyle(
      seqs(['chr1', 'chr2', 'chrX', 'chr1_KI270706v1_random', 'chrM', 'chrEBV']),
    );
    expect(naming).toEqual({ style: 'ucsc', examples: ['chr1', 'chrM'] });
    expect(namingLabel(naming!)).toBe('UCSC-style names (chr1, chrM)');
  });

  it('reads Ensembl names', () => {
    const naming = namingStyle(seqs(['1', '2', 'X', 'Y', 'MT', 'GL456210.1']));
    expect(namingLabel(naming!)).toBe('Ensembl-style names (1, MT)');
  });

  it('reads accession names', () => {
    expect(namingStyle(seqs(['NC_000001.11', 'NC_000002.12', 'NT_187361.1']))?.style).toBe(
      'accession',
    );
  });

  it('decides by the longest sequences, not the many scaffolds', () => {
    const main = seqs(['chr1', 'chr2', 'chr3']);
    const scaffolds = Array.from({ length: 50 }, (_, i) => ({
      name: `KI2707${String(i).padStart(2, '0')}.1`,
      length: 100 + i,
    }));
    expect(namingStyle([...scaffolds, ...main])?.style).toBe('ucsc');
  });

  it('gives nothing for names in no known style, or none at all', () => {
    expect(namingStyle(seqs(['SUPER_1', 'SUPER_2', 'scaffold_3']))).toBeNull();
    expect(namingStyle([])).toBeNull();
  });
});
