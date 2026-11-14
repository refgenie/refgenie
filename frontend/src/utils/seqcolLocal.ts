/**
 * A plain-language summary of how a user's FASTA file relates to a genome,
 * from the two lists of sequences.
 */

/** One sequence of a collection, in collection order. */
export interface SeqcolRow {
  name: string;
  length: number;
  /** Bare GA4GH digest, no `SQ.` prefix. */
  sha512t24u: string;
}

export interface FileSummary {
  /** One short line: how the file relates to the genome. */
  headline: string;
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

  const coords = (rows: readonly SeqcolRow[]) =>
    new Set(rows.map((row) => `${row.name}\t${row.length}`));
  const genomeCoords = coords(genome);
  const fileCoords = coords(file);
  if (
    genomeCoords.size === fileCoords.size &&
    [...fileCoords].every((pair) => genomeCoords.has(pair))
  ) {
    details.push(
      'The names and lengths match, so positions (in a BED file, say) mean the same in both.',
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
  return { headline, details };
}
