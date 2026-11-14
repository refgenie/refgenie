import { describe, expect, it } from 'vitest';
import { namingStyle } from './namingStyle';

const seqs = (names: string[]) =>
  names.map((name, i) => ({ name, length: 1_000_000 - i * 1000 }));

describe('namingStyle', () => {
  it('reads the longest sequence name', () => {
    expect(namingStyle(seqs(['chr1', 'chr2', 'chrM']))).toBe('UCSC-style names (chr1)');
    expect(namingStyle(seqs(['1', '2', 'MT']))).toBe('Ensembl-style names (1)');
    expect(namingStyle([{ name: 'KI270706.1', length: 9 }, ...seqs(['X'])])).toBe(
      'Ensembl-style names (X)',
    );
  });

  it('gives nothing for other names, or none at all', () => {
    expect(namingStyle(seqs(['NC_000001.11', 'chr2']))).toBeNull();
    expect(namingStyle([])).toBeNull();
  });
});
