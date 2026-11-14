import { describe, expect, it } from 'vitest';
import {
  compareSeqcols,
  comparisonArrays,
  summarizeFileComparison,
  type SeqcolRow,
} from './seqcolLocal';
import { interpretCoordinateSystem } from './seqcolCompare';

const row = (name: string, length: number, seq: string): SeqcolRow => ({
  name,
  length,
  sha512t24u: seq,
});

// refget's demo_fasta/demo0.fa and demo1.fa.gz; the expected counts are
// refget's demo_fasta/compare-0vs1.json, from the reference implementation.
const DEMO0 = [row('chrX', 8, 'x0'), row('chr1', 4, 'a0'), row('chr2', 4, 'b0')];
const DEMO1 = [row('chr1', 4, 'a1'), row('chr2', 4, 'b1')];

const GENOME = [row('chr1', 100, 's1'), row('chr2', 80, 's2'), row('chrM', 16, 'sM')];

describe('compareSeqcols', () => {
  it('matches the reference comparison of demo0 with demo1', () => {
    const c = compareSeqcols(comparisonArrays(DEMO0), comparisonArrays(DEMO1), {
      a: 'A',
      b: 'B',
    });
    const e = c.array_elements;
    expect(e.a_count).toMatchObject({ names: 3, lengths: 3, sequences: 3 });
    expect(e.b_count).toMatchObject({ names: 2, lengths: 2, sequences: 2 });
    expect(e.a_and_b_count).toMatchObject({ names: 2, lengths: 2, sequences: 0 });
    expect(e.a_and_b_same_order).toMatchObject({
      names: true,
      lengths: true,
      sequences: null,
    });
    expect(c.attributes).toEqual({
      a_only: [],
      b_only: [],
      a_and_b: [
        'lengths',
        'name_length_pairs',
        'names',
        'sequences',
        'sorted_sequences',
      ],
    });
  });

  it('reports the same five attributes as the server, and identity', () => {
    const arrays = comparisonArrays(GENOME);
    const c = compareSeqcols(arrays, arrays, { a: 'G', b: 'G' });
    expect(interpretCoordinateSystem(c).level).toBe('same');
    for (const attr of c.attributes.a_and_b) {
      expect(c.array_elements.a_and_b_count[attr]).toBe(3);
      expect(c.array_elements.a_and_b_same_order[attr]).toBe(true);
    }
  });

  it('sees a reorder in sequences but not in sorted_sequences', () => {
    const file = [GENOME[1], GENOME[0], GENOME[2]];
    const c = compareSeqcols(comparisonArrays(GENOME), comparisonArrays(file), {
      a: 'G',
      b: 'F',
    });
    expect(c.array_elements.a_and_b_same_order.sequences).toBe(false);
    expect(c.array_elements.a_and_b_same_order.sorted_sequences).toBe(true);
  });

  it('leaves order unknown with unbalanced duplicates', () => {
    const c = compareSeqcols(
      { x: ['a', 'b', 'a'] },
      { x: ['a', 'b'] },
      { a: 'A', b: 'B' },
    );
    expect(c.array_elements.a_and_b_count.x).toBe(2);
    expect(c.array_elements.a_and_b_same_order.x).toBeNull();
  });

  it('lists attributes only one side has', () => {
    const c = compareSeqcols(
      { x: ['1'], y: ['2'] },
      { x: ['1'], z: ['3'] },
      { a: 'A', b: 'B' },
    );
    expect(c.attributes).toEqual({ a_only: ['y'], b_only: ['z'], a_and_b: ['x'] });
  });
});

describe('summarizeFileComparison', () => {
  const digests = { genome: 'G', file: 'F' };

  it('says when the file is this genome', () => {
    const s = summarizeFileComparison(GENOME, GENOME, { genome: 'G', file: 'G' });
    expect(s.identical).toBe(true);
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
