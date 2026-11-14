import { describe, expect, it } from 'vitest';
import { alphabetCounts } from './alphabets';
import type { SequenceRow } from './types';

const rows = (...specs: [string | null, number][]) =>
  specs.map(([alphabet, length]) => ({ alphabet, length }) as unknown as SequenceRow);

describe('alphabetCounts', () => {
  it('counts sequences and bases for each alphabet, most sequences first', () => {
    expect(
      alphabetCounts(rows(['dna2bit', 10], ['dna3bit', 1000], ['dna2bit', 20], [null, 5])),
    ).toEqual([
      { alphabet: 'dna2bit', count: 2, bases: 30 },
      { alphabet: 'dna3bit', count: 1, bases: 1000 },
      { alphabet: 'unknown', count: 1, bases: 5 },
    ]);
  });

  it('folds a fifth alphabet or more into other, summing both measures', () => {
    const counts = alphabetCounts(
      rows(['a', 1], ['a', 1], ['a', 1], ['b', 1], ['b', 1], ['c', 1], ['d', 4], ['e', 6]),
    );
    expect(counts.map((c) => c.alphabet)).toEqual(['a', 'b', 'c', 'other']);
    expect(counts[3]).toEqual({ alphabet: 'other', count: 2, bases: 10 });
  });
});
