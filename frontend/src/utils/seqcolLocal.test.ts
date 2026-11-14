import { describe, expect, it } from 'vitest';
import { summarizeFileComparison, type SeqcolRow } from './seqcolLocal';

const row = (name: string, length: number, seq: string): SeqcolRow => ({
  name,
  length,
  sha512t24u: seq,
});

const GENOME = [row('chr1', 100, 's1'), row('chr2', 80, 's2'), row('chrM', 16, 'sM')];

describe('summarizeFileComparison', () => {
  const digests = { genome: 'G', file: 'F' };

  it('says when the file is this genome', () => {
    const s = summarizeFileComparison(GENOME, GENOME, { genome: 'G', file: 'G' });
    expect(s.headline).toBe('Your file is this genome');
  });

  it('same sequences, different names', () => {
    const file = [row('1', 100, 's1'), row('2', 80, 's2'), row('MT', 16, 'sM')];
    const s = summarizeFileComparison(GENOME, file, digests);
    expect(s.headline).toBe('Same sequences, different names');
    expect(s.details).toContain(
      '3 shared sequences have a different name in your file (chr1 is 1, chr2 is 2, chrM is MT).',
    );
  });

  it('a renamed and a dropped sequence', () => {
    const file = [row('chr1', 100, 's1'), row('MT', 16, 'sM')];
    const s = summarizeFileComparison(GENOME, file, digests);
    expect(s.headline).toBe('Your file is part of this genome');
    expect(s.details).toEqual([
      '2 sequences in your file are also in this genome.',
      'Your file is missing 1 sequence of this genome (chr2).',
      '1 shared sequence has a different name in your file (chrM is MT).',
    ]);
  });

  it('extra sequences, and a name reused for other bases', () => {
    const file = [
      ...GENOME.slice(0, 2),
      row('chrM', 16, 'other'),
      row('chrX', 9, 'x'),
      row('chrY', 9, 'y'),
      row('chrZ', 9, 'z'),
      row('chrW', 9, 'w'),
    ];
    const s = summarizeFileComparison(GENOME, file, digests);
    expect(s.headline).toBe('Some sequences in common');
    expect(s.details).toContain(
      'Your file has 5 sequences that are not in this genome (for example chrM, chrX, chrY).',
    );
    expect(s.details).toContain(
      '1 name is used in both, but for different sequences (chrM).',
    );
  });

  it('says when the names and lengths match though the bases differ', () => {
    const file = [row('chr1', 100, 's1'), row('chr2', 80, 'patched'), row('chrM', 16, 'sM')];
    expect(summarizeFileComparison(GENOME, file, digests).details).toContain(
      'The names and lengths match, so positions (in a BED file, say) mean the same in both.',
    );
  });

  it('a superset, a reorder, and nothing in common', () => {
    expect(
      summarizeFileComparison(GENOME, [...GENOME, row('chrX', 9, 'x')], digests)
        .headline,
    ).toBe('Your file has this whole genome, and more');
    expect(
      summarizeFileComparison(GENOME, [GENOME[2], GENOME[0], GENOME[1]], digests)
        .headline,
    ).toBe('Same sequences and names, in a different order');
    expect(summarizeFileComparison(GENOME, [row('x', 1, 'x')], digests).headline).toBe(
      'No sequences in common',
    );
  });
});
