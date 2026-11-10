// Ported from repos/refget/frontend/src/pages/StoreCollection.tsx (summary stats) @ c0d3e09; keep the two in step.
// The sequence count is on the Sequences tab, so this shows the total size and the alphabets.

import { alphabetCounts } from './alphabets';
import type { AlphabetPart } from './alphabets';
import { formatBp } from '../../utils/format';
import type { SequenceRow } from './types';

export interface SequenceStatsProps {
  sequences: SequenceRow[];
}

/** Each part's left edge and width in percent, with a small gap between parts. */
function segments(parts: AlphabetPart[], measure: 'count' | 'bases', total: number) {
  const gap = 0.6;
  const usable = 100 - gap * (parts.length - 1);
  let x = 0;
  return parts.map((part) => {
    const width = Math.max((usable * part[measure]) / total, 0.5);
    const segment = { alphabet: part.alphabet, x, width };
    x += width + gap;
    return segment;
  });
}

const plural = (n: number) =>
  `${n.toLocaleString()} ${n === 1 ? 'sequence' : 'sequences'}`;

/** Total size and alphabets, compact, beside the filter box. */
export function SequenceStats({ sequences }: SequenceStatsProps) {
  const bases = sequences.reduce((sum, s) => sum + s.length, 0);
  return (
    <div className="rg-seq-summary">
      <p className="rg-seq-summary__item">
        <span className="rg-seq-summary__label">Total size</span>
        <span title={`${bases.toLocaleString()} bp`}>{formatBp(bases)}</span>
      </p>
      <Alphabets sequences={sequences} />
    </div>
  );
}

const pct = (part: number, whole: number) => {
  const p = (100 * part) / whole;
  return p > 0 && p < 1 ? '<1%' : `${Math.round(p)}%`;
};

/**
 * How the genome splits across alphabets, two ways: by sequences and by
 * bases. They often differ a lot: long chromosomes with N runs need dna3bit,
 * while many short contigs fit in dna2bit. One alphabet is a word.
 */
function Alphabets({ sequences }: SequenceStatsProps) {
  const parts = alphabetCounts(sequences);
  if (parts.length === 0) return null;
  if (parts.length === 1) {
    return (
      <p className="rg-seq-summary__item">
        <span className="rg-seq-summary__label">Alphabet</span>
        <code>{parts[0].alphabet}</code>
      </p>
    );
  }
  const total = sequences.length;
  const bases = parts.reduce((sum, p) => sum + p.bases, 0);
  const bars: { measure: 'count' | 'bases'; label: string; whole: number }[] = [
    { measure: 'count', label: 'sequences', whole: total },
    { measure: 'bases', label: 'bases', whole: bases },
  ];
  return (
    <figure className="rg-alphabets rg-seq-summary__item">
      <figcaption className="rg-seq-summary__label">Alphabets</figcaption>
      <div className="rg-alphabets__body">
        {bars.map(({ measure, label, whole }) => (
          <div key={measure} className="rg-alphabets__row">
            <span className="rg-alphabets__measure">{label}</span>
            <svg className="rg-alphabets__bar" aria-hidden="true">
              {segments(parts, measure, whole).map(({ alphabet, x, width }, i) => (
                <rect
                  key={alphabet}
                  x={`${x}%`}
                  width={`${width}%`}
                  height="100%"
                  className={`rg-alphabets__seg rg-series-${i + 1}`}
                >
                  <title>
                    {measure === 'count'
                      ? `${alphabet}: ${plural(parts[i].count)} (${pct(parts[i].count, whole)})`
                      : `${alphabet}: ${formatBp(parts[i].bases)} (${pct(parts[i].bases, whole)})`}
                  </title>
                </rect>
              ))}
            </svg>
          </div>
        ))}
        <ul className="rg-alphabets__list">
          {parts.map((part, i) => (
            <li
              key={part.alphabet}
              className="rg-alphabets__item"
              title={`${plural(part.count)}, ${part.bases.toLocaleString()} bp`}
            >
              <span className={`rg-alphabets__swatch rg-series-${i + 1}`} aria-hidden="true" />
              <code>{part.alphabet}</code>
              <span className="rg-muted">
                {pct(part.count, total)} of seqs · {pct(part.bases, bases)} of bp
              </span>
            </li>
          ))}
        </ul>
      </div>
    </figure>
  );
}
