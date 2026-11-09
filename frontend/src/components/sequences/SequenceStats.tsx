// Ported from repos/refget/frontend/src/pages/StoreCollection.tsx (summary stats) @ c0d3e09; keep the two in step.

import type { SequenceRow } from './types';

export interface SequenceStatsProps {
  sequences: SequenceRow[];
}

/** One line: sequence count, total bases and how many sequences use each alphabet. */
export function SequenceStats({ sequences }: SequenceStatsProps) {
  const totalBases = sequences.reduce((sum, s) => sum + s.length, 0);
  const alphabetCounts = new Map<string, number>();
  for (const s of sequences) {
    const alphabet = s.alphabet ?? 'unknown';
    alphabetCounts.set(alphabet, (alphabetCounts.get(alphabet) ?? 0) + 1);
  }
  const alphabets = [...alphabetCounts]
    .map(([alphabet, count]) => `${alphabet} ${count.toLocaleString()}`)
    .join(', ');

  return (
    <p className="rg-stat-line">
      <span>
        {sequences.length.toLocaleString()}{' '}
        {sequences.length === 1 ? 'sequence' : 'sequences'}
      </span>
      <span aria-hidden="true"> · </span>
      <span>{totalBases.toLocaleString()} bp</span>
      {alphabets && (
        <>
          <span aria-hidden="true"> · </span>
          <span>alphabets: {alphabets}</span>
        </>
      )}
    </p>
  );
}

/** A placeholder the same height as {@link SequenceStats}. */
export function SequenceStatsSkeleton() {
  return (
    <p className="rg-stat-line" aria-hidden="true">
      <span className="rg-skeleton rg-skeleton--inline rg-skeleton--w24" />
    </p>
  );
}
