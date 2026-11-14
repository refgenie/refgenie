import { describe, expect, it } from 'vitest';
import { interpretAttribute, interpretCoordinateSystem } from './seqcolCompare';
import type { SeqcolComparison } from '../types/api';

const ATTRS = ['names', 'lengths', 'name_length_pairs', 'sequences'] as const;

/** Build a comparison from [a, b, shared] counts per attribute. */
function cmp(
  counts: Record<(typeof ATTRS)[number], [number, number, number]>,
  sameOrder: boolean | null = true,
): SeqcolComparison {
  const side = (i: 0 | 1 | 2) =>
    Object.fromEntries(ATTRS.map((attr) => [attr, counts[attr][i]]));
  return {
    digests: { a: 'A', b: 'B' },
    attributes: { a_only: [], b_only: [], a_and_b: [...ATTRS] },
    array_elements: {
      a_count: side(0),
      b_count: side(1),
      a_and_b_count: side(2),
      a_and_b_same_order: Object.fromEntries(ATTRS.map((attr) => [attr, sameOrder])),
    },
  };
}

const all = (a: number, b: number, shared: number) => ({
  names: [a, b, shared] as [number, number, number],
  lengths: [a, b, shared] as [number, number, number],
  name_length_pairs: [a, b, shared] as [number, number, number],
  sequences: [a, b, shared] as [number, number, number],
});

const EMOJI = /\p{Extended_Pictographic}/u;

function allMessages(c: SeqcolComparison): string[] {
  return [...interpretCoordinateSystem(c).messages, ...ATTRS.map((a) => interpretAttribute(c, a))];
}

describe('seqcol comparison interpretation', () => {
  it('identical', () => {
    const c = cmp(all(3, 3, 3));
    expect(interpretCoordinateSystem(c).level).toBe('same');
    expect(interpretAttribute(c, 'names')).toMatch(/identical.*same order/);
  });

  it('subset', () => {
    const c = cmp(all(1, 2, 1));
    const verdict = interpretCoordinateSystem(c);
    expect(verdict.level).toBe('subset');
    expect(verdict.messages[0]).toMatch(/This genome's coordinate system is a subset/);
    expect(interpretAttribute(c, 'sequences')).toMatch(/all 1 sequences of this genome, and 1 more/);
  });

  it('names the other side as asked', () => {
    const c = cmp(all(2, 1, 1));
    expect(interpretCoordinateSystem(c, 'your file').messages[0]).toBe(
      "Your file's coordinate system is a subset of this genome's.",
    );
    expect(interpretAttribute(c, 'names', 'your file')).toBe(
      'This genome has all 1 names of your file, and 1 more.',
    );
  });

  it('partial overlap', () => {
    const c = cmp(all(3, 3, 1));
    expect(interpretCoordinateSystem(c).level).toBe('overlap');
    expect(interpretAttribute(c, 'names')).toMatch(/partly overlap/);
  });

  it('disjoint', () => {
    const c = cmp(all(2, 2, 0));
    expect(interpretCoordinateSystem(c).level).toBe('disjoint');
    expect(interpretAttribute(c, 'lengths')).toMatch(/disjoint/);
  });

  it('name swap', () => {
    const c = cmp({
      names: [2, 2, 2],
      lengths: [2, 2, 2],
      name_length_pairs: [2, 2, 0],
      sequences: [2, 2, 2],
    });
    const verdict = interpretCoordinateSystem(c);
    expect(verdict.level).toBe('swap');
    expect(verdict.messages.join(' ')).toMatch(/name swap/);
  });

  it('never uses emoji', () => {
    for (const c of [cmp(all(3, 3, 3), false), cmp(all(1, 2, 1)), cmp(all(2, 2, 0))]) {
      for (const message of allMessages(c)) expect(message).not.toMatch(EMOJI);
    }
  });
});
