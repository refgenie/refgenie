/**
 * The GA4GH seqcol comparison, computed in the browser from two lists of
 * sequences, plus a plain-language summary of how a user's FASTA file relates
 * to a genome.
 *
 * `compareSeqcols` follows the seqcol spec's comparison function as gtars
 * implements it (`gtars-refget` `compare_arrays` / `compare_elements`), on the
 * same five attributes the server's `/compare` endpoint reports. The published
 * gtars wasm has no comparison export, hence this port.
 */

import type { SeqcolComparison } from '../types/api';

/** One sequence of a collection, in collection order. */
export interface SeqcolRow {
  name: string;
  length: number;
  /** Bare GA4GH digest, no `SQ.` prefix. */
  sha512t24u: string;
}

export type SeqcolArrays = Record<string, string[]>;

/**
 * The level-2 arrays a comparison reads. `name_length_pairs` elements are
 * `{"length":N,"name":"..."}` JSON; both sides are built here, so the two
 * encodings always agree.
 */
export function comparisonArrays(rows: readonly SeqcolRow[]): SeqcolArrays {
  const sequences = rows.map((row) => `SQ.${row.sha512t24u}`);
  return {
    names: rows.map((row) => row.name),
    lengths: rows.map((row) => String(row.length)),
    sequences,
    sorted_sequences: [...sequences].sort(),
    name_length_pairs: rows.map((row) =>
      JSON.stringify({ length: row.length, name: row.name }),
    ),
  };
}

/**
 * Elements of `a` found in `b` and of `b` found in `a`: the overlap is the
 * smaller of the two filtered lengths, and the order agrees when both
 * filtered lists are equal. Order is `null` with fewer than two shared
 * elements or unbalanced duplicates.
 */
function compareElements(a: string[], b: string[]): [number, boolean | null] {
  const inA = new Set(a);
  const inB = new Set(b);
  const fa = a.filter((x) => inB.has(x));
  const fb = b.filter((x) => inA.has(x));
  const overlap = Math.min(fa.length, fb.length);
  if (overlap < 2 || fa.length !== fb.length) return [overlap, null];
  return [overlap, fa.every((x, i) => x === fb[i])];
}

/** The seqcol comparison of collection `a` with collection `b`. */
export function compareSeqcols(
  a: SeqcolArrays,
  b: SeqcolArrays,
  digests: { a: string; b: string },
): SeqcolComparison {
  const keys = [...new Set([...Object.keys(a), ...Object.keys(b)])].sort();
  const result: SeqcolComparison = {
    digests,
    attributes: { a_only: [], b_only: [], a_and_b: [] },
    array_elements: {
      a_count: {},
      b_count: {},
      a_and_b_count: {},
      a_and_b_same_order: {},
    },
  };
  for (const key of keys) {
    const inA = key in a;
    const inB = key in b;
    if (inA) result.array_elements.a_count[key] = a[key].length;
    if (inB) result.array_elements.b_count[key] = b[key].length;
    if (inA && inB) {
      result.attributes.a_and_b.push(key);
      const [overlap, order] = compareElements(a[key], b[key]);
      result.array_elements.a_and_b_count[key] = overlap;
      result.array_elements.a_and_b_same_order[key] = order;
    } else if (inA) {
      result.attributes.a_only.push(key);
    } else {
      result.attributes.b_only.push(key);
    }
  }
  return result;
}

export interface FileSummary {
  /** One short line: how the file relates to the genome. */
  headline: string;
  identical: boolean;
  /** Plain sentences with up to a few examples each. */
  details: string[];
}

const EXAMPLES = 3;

function plural(n: number, one: string, many = `${one}s`): string {
  return `${n.toLocaleString()} ${n === 1 ? one : many}`;
}

function examples(items: string[]): string {
  if (items.length === 0) return '';
  const shown = items.slice(0, EXAMPLES).join(', ');
  return items.length > EXAMPLES ? ` (for example ${shown})` : ` (${shown})`;
}

function firstBy<K>(rows: readonly SeqcolRow[], key: (row: SeqcolRow) => K) {
  const map = new Map<K, SeqcolRow>();
  for (const row of rows) if (!map.has(key(row))) map.set(key(row), row);
  return map;
}

/**
 * Describe a user's FASTA file (`file`) against this genome (`genome`) in
 * words: which sequences are shared, extra, missing, renamed, or reordered.
 * Sequences match by content digest; names only label them.
 */
export function summarizeFileComparison(
  genome: readonly SeqcolRow[],
  file: readonly SeqcolRow[],
  digests: { genome: string; file: string },
): FileSummary {
  if (digests.genome === digests.file) {
    return {
      headline: 'Your file is this genome',
      identical: true,
      details: ['Same sequences, same names, same order. The digests match.'],
    };
  }
  const genomeBySeq = firstBy(genome, (row) => row.sha512t24u);
  const fileBySeq = firstBy(file, (row) => row.sha512t24u);
  const genomeByName = firstBy(genome, (row) => row.name);

  const extra = [...fileBySeq.values()].filter(
    (row) => !genomeBySeq.has(row.sha512t24u),
  );
  const missing = [...genomeBySeq.values()].filter(
    (row) => !fileBySeq.has(row.sha512t24u),
  );
  const shared = [...fileBySeq.values()].filter((row) =>
    genomeBySeq.has(row.sha512t24u),
  );
  const renamed = shared.filter(
    (row) => genomeBySeq.get(row.sha512t24u)?.name !== row.name,
  );
  const changed = extra.filter((row) => genomeByName.has(row.name));

  const details: string[] = [];
  if (shared.length > 0) {
    details.push(
      `${plural(shared.length, 'sequence')} in your file ${shared.length === 1 ? 'is' : 'are'} also in this genome.`,
    );
  }
  if (extra.length > 0) {
    details.push(
      `Your file has ${plural(extra.length, 'sequence')} that ${extra.length === 1 ? 'is' : 'are'} not in this genome${examples(extra.map((row) => row.name))}.`,
    );
  }
  if (missing.length > 0) {
    details.push(
      `Your file is missing ${plural(missing.length, 'sequence')} of this genome${examples(missing.map((row) => row.name))}.`,
    );
  }
  if (renamed.length > 0) {
    details.push(
      `${plural(renamed.length, 'shared sequence')} ${renamed.length === 1 ? 'has' : 'have'} a different name in your file${examples(
        renamed.map((row) => `${genomeBySeq.get(row.sha512t24u)?.name} is ${row.name}`),
      )}.`,
    );
  }
  if (changed.length > 0) {
    details.push(
      `${plural(changed.length, 'name')} ${changed.length === 1 ? 'is' : 'are'} used in both, but for different sequences${examples(changed.map((row) => row.name))}.`,
    );
  }

  const sameSet = extra.length === 0 && missing.length === 0;
  let headline: string;
  if (sameSet && renamed.length > 0) headline = 'Same sequences, different names';
  else if (sameSet) {
    headline =
      genome.length === file.length
        ? 'Same sequences and names, in a different order'
        : 'Same sequences, some repeated';
  } else if (shared.length === 0) headline = 'No sequences in common';
  else if (extra.length === 0) headline = 'Your file is part of this genome';
  else if (missing.length === 0) headline = 'Your file has this whole genome, and more';
  else headline = 'Some sequences in common';
  return { headline, identical: false, details };
}
