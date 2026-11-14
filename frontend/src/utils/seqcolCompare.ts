// Ported from repos/refget/frontend/src/pages/ComparisonView.tsx @ c0d3e09; keep the two in step.
// Plain-language reading of a GA4GH seqcol comparison. Collection A is "this
// genome", collection B "the other genome" (or `other`, e.g. "your file"). No
// emoji: a `level` picks the tone.

import type { SeqcolComparison } from '../types/api';

export type CompareLevel = 'same' | 'subset' | 'overlap' | 'disjoint' | 'swap';

const OTHER = 'the other genome';

function cap(text: string): string {
  return text.charAt(0).toUpperCase() + text.slice(1);
}

function count(c: SeqcolComparison, side: keyof SeqcolComparison['array_elements'], attr: string) {
  const value = c.array_elements[side]?.[attr];
  return typeof value === 'number' ? value : undefined;
}

/** How the two coordinate systems (name and length pairs) relate. */
export function interpretCoordinateSystem(
  c: SeqcolComparison,
  other = OTHER,
): {
  level: CompareLevel;
  messages: string[];
} {
  const nlpA = count(c, 'a_count', 'name_length_pairs');
  const nlpB = count(c, 'b_count', 'name_length_pairs');
  const nlpBoth = count(c, 'a_and_b_count', 'name_length_pairs');
  if (nlpA === undefined || nlpB === undefined || nlpBoth === undefined) {
    return { level: 'disjoint', messages: ['Unable to interpret: incomplete comparison data.'] };
  }
  const missing = (attr: string, side: 'a_count' | 'b_count') =>
    (count(c, side, attr) ?? 0) - (count(c, 'a_and_b_count', attr) ?? 0);
  const nlpANotB = nlpA - nlpBoth;
  const nlpBNotA = nlpB - nlpBoth;

  if (nlpANotB === 0 && nlpBNotA === 0) {
    return {
      level: 'same',
      messages: ['The coordinate systems are identical.', 'Names match. Lengths match.'],
    };
  }
  if (nlpANotB === 0) {
    return {
      level: 'subset',
      messages: [`This genome's coordinate system is a subset of ${other}'s.`],
    };
  }
  if (nlpBNotA === 0) {
    return {
      level: 'subset',
      messages: [`${cap(other)}'s coordinate system is a subset of this genome's.`],
    };
  }
  if (nlpBoth !== 0) {
    return { level: 'overlap', messages: ['The coordinate systems partly overlap.'] };
  }

  const lengthsMatch = missing('lengths', 'a_count') === 0 && missing('lengths', 'b_count') === 0;
  const namesMatch = missing('names', 'a_count') === 0 && missing('names', 'b_count') === 0;
  const messages = ['The coordinate systems are disjoint.'];
  if (lengthsMatch && namesMatch) {
    messages.push('The same names and lengths are paired differently (a name swap).');
    return { level: 'swap', messages };
  }
  if (lengthsMatch) messages.push("Lengths match. Names don't match.");
  else if (namesMatch) messages.push("Names match. Lengths don't match.");
  return { level: 'disjoint', messages };
}

/** One sentence (or two) on how one attribute's contents compare. */
export function interpretAttribute(
  c: SeqcolComparison,
  attribute: string,
  other = OTHER,
): string {
  const nA = count(c, 'a_count', attribute);
  const nB = count(c, 'b_count', attribute);
  const both = count(c, 'a_and_b_count', attribute);
  if (nA === undefined || nB === undefined || both === undefined) {
    return `Unable to interpret: missing ${attribute} data.`;
  }
  const order = c.array_elements.a_and_b_same_order?.[attribute];
  const messages: string[] = [];

  if (both === nA && both === nB) {
    messages.push(`The ${attribute} contents are identical.`);
    if (order === true) messages.push('The elements are in the same order.');
    else if (order === false) messages.push('The elements are in a different order.');
  }
  if (both === nA && both < nB) {
    messages.push(
      `${cap(other)} has all ${nA} ${attribute} of this genome, and ${nB - both} more.`,
    );
  }
  if (both === nB && both < nA) {
    messages.push(`This genome has all ${nB} ${attribute} of ${other}, and ${nA - both} more.`);
  }
  if (both === 0 && (nA > 0 || nB > 0)) {
    messages.push(`The ${attribute} contents are disjoint.`);
  } else if (both < nA && both < nB) {
    messages.push(
      `The ${attribute} contents partly overlap: some are shared, and some are unique to each.`,
    );
  }
  return messages.join(' ');
}
