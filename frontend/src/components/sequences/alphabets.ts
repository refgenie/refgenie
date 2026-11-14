import type { SequenceRow } from './types';

/** Categorical slots in fixed order; a fifth alphabet or more folds into "other". */
const MAX_SLOTS = 4;

export interface AlphabetPart {
  alphabet: string;
  /** Sequences that use it. */
  count: number;
  /** Total bases in those sequences. */
  bases: number;
}

/**
 * Alphabets with their sequence counts and total bases, most sequences first,
 * the tail folded into "other". The order (and so each alphabet's colour) is
 * the same for both measures.
 */
export function alphabetCounts(sequences: SequenceRow[]): AlphabetPart[] {
  const parts = new Map<string, AlphabetPart>();
  for (const s of sequences) {
    const alphabet = s.alphabet ?? 'unknown';
    const part = parts.get(alphabet) ?? { alphabet, count: 0, bases: 0 };
    part.count += 1;
    part.bases += s.length;
    parts.set(alphabet, part);
  }
  const sorted = [...parts.values()].sort(
    (a, b) => b.count - a.count || a.alphabet.localeCompare(b.alphabet),
  );
  if (sorted.length <= MAX_SLOTS) return sorted;
  const head = sorted.slice(0, MAX_SLOTS - 1);
  const tail = sorted.slice(MAX_SLOTS - 1);
  return [
    ...head,
    {
      alphabet: 'other',
      count: tail.reduce((sum, a) => sum + a.count, 0),
      bases: tail.reduce((sum, a) => sum + a.bases, 0),
    },
  ];
}
