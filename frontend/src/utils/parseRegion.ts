/**
 * Parse a region typed or pasted by a user, samtools/UCSC style:
 *
 *   "chr1:1,000,000-1,000,500"  bases 1,000,000 through 1,000,500 (501 bases)
 *   "chr1:1000000-1000500"      the same
 *   "chr1:1000000"              the single base at 1,000,000
 *   "chr1:1000000-"             from 1,000,000 to the end of chr1
 *   "chrM"                      the whole sequence
 *
 * Input coordinates are 1-based and inclusive, like `samtools faidx`. The
 * result is 0-based and half-open, the form the store reads take.
 *
 * Adapted from the RefgetStore Explorer's parseRegion (which reads 0-based,
 * half-open coordinates and assumes the sequence is already chosen). Here the
 * sequence name is part of the input and is matched against sequence names
 * and their aliases.
 */

/** The most bases one in-browser extraction may return. */
export const MAX_EXTRACT_BASES = 100_000;

export interface RegionSequence {
  name: string;
  length: number;
  aliases?: { alias: string }[];
}

export interface ParsedRegion {
  /** The sequence's own name (an alias in the input is resolved to it). */
  name: string;
  /** 0-based, inclusive. */
  start: number;
  /** 0-based, exclusive. */
  end: number;
  /** True if the end was past the sequence and was cut back to its length. */
  clamped: boolean;
  /** The region in 1-based, inclusive form, for a FASTA header or a link. */
  label: string;
}

export class RegionError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'RegionError';
  }
}

const fmt = (n: number) => n.toLocaleString('en-US');

/**
 * A region that exists in this genome, for placeholders, error hints and the
 * "Try an example" button: 501 bases from the longest sequence, at 1,000,000
 * when it is long enough, else from its middle.
 */
export function exampleRegion(sequences: readonly RegionSequence[]): string | null {
  let longest: RegionSequence | null = null;
  for (const seq of sequences) if (!longest || seq.length > longest.length) longest = seq;
  if (!longest || longest.length < 1) return null;
  const len = longest.length;
  const start = len > 1_000_500 ? 1_000_000 : Math.max(1, Math.floor(len / 2) - 250);
  const end = Math.min(len, start + 500);
  return `${longest.name}:${fmt(start)}-${fmt(end)}`;
}
/**
 * `name:start[-[end]]`. Spaces may sit around ':' and '-' and between digit
 * groups, but the name keeps its own spaces.
 */
const NUM = '[\\d,_]+(?:\\s+[\\d,_]+)*';
const COORDS = new RegExp(`^(.*\\S)\\s*:\\s*(${NUM})(\\s*-\\s*(${NUM}|))?$`);

function toNumber(raw: string): number {
  return Number(raw.replace(/[,_\s]/g, ''));
}

/** Find a sequence by exact name, then by exact alias. */
export function findSequence<T extends RegionSequence>(
  name: string,
  sequences: readonly T[],
): T | null {
  return (
    sequences.find((seq) => seq.name === name) ??
    sequences.find((seq) => seq.aliases?.some((a) => a.alias === name)) ??
    null
  );
}

/**
 * Up to three sequence names close to `name`: the same name in another case,
 * with or without a "chr" prefix, or M for MT and back.
 */
export function suggestSequences(
  name: string,
  sequences: readonly RegionSequence[],
): string[] {
  const bare = name.toLowerCase().replace(/^chr/, '');
  const wanted = new Set([bare, `chr${bare}`]);
  if (bare === 'm' || bare === 'mt') {
    for (const v of ['m', 'mt', 'chrm', 'chrmt']) wanted.add(v);
  }
  const out: string[] = [];
  for (const seq of sequences) {
    const keys = [seq.name, ...(seq.aliases ?? []).map((a) => a.alias)];
    if (keys.some((key) => wanted.has(key.toLowerCase())) && !out.includes(seq.name)) {
      out.push(seq.name);
      if (out.length === 3) break;
    }
  }
  return out;
}

function unknownName(name: string, sequences: readonly RegionSequence[]): RegionError {
  const close = suggestSequences(name, sequences);
  const hint = close.length
    ? ` Did you mean ${close.map((n) => `"${n}"`).join(' or ')}?`
    : '';
  return new RegionError(`This genome has no sequence named "${name}".${hint}`);
}

export function parseRegion(
  text: string,
  sequences: readonly RegionSequence[],
  maxBases = MAX_EXTRACT_BASES,
): ParsedRegion {
  const input = String(text ?? '').trim();
  const example = exampleRegion(sequences);
  const forExample = example ? `, for example ${example}` : '';
  if (!input) throw new RegionError(`Enter a region${forExample}.`);

  // A whole name first, so names that contain ":" (e.g. HLA-A*01:01:01:01)
  // are not misread as coordinates.
  let seq = findSequence(input, sequences);
  let startRaw: number | null = null;
  let endRaw: number | null = null;
  if (!seq) {
    const match = COORDS.exec(input);
    if (!match) {
      if (/[:]/.test(input)) {
        throw new RegionError(
          `Cannot read "${input}" as a region.${example ? ` Try ${example}.` : ''}`,
        );
      }
      throw unknownName(input, sequences);
    }
    const name = match[1];
    seq = findSequence(name, sequences);
    if (!seq) throw unknownName(name, sequences);
    startRaw = toNumber(match[2]);
    if (match[3] === undefined) endRaw = startRaw;
    else if (match[4] !== '') endRaw = toNumber(match[4]);
  }

  const length = seq.length;
  const startOne = startRaw ?? 1;
  const endOne = endRaw ?? length;
  if (!Number.isSafeInteger(startOne) || !Number.isSafeInteger(endOne)) {
    throw new RegionError(
          `Cannot read "${input}" as a region.${example ? ` Try ${example}.` : ''}`,
        );
  }
  if (startOne < 1) {
    throw new RegionError('Positions start at 1 (1-based, like samtools faidx).');
  }
  if (endOne < startOne) {
    throw new RegionError(
      `The end (${endOne.toLocaleString()}) is before the start (${startOne.toLocaleString()}).`,
    );
  }
  if (startOne > length) {
    throw new RegionError(
      `${seq.name} has ${length.toLocaleString()} bases, so it has no position ` +
        `${startOne.toLocaleString()}.`,
    );
  }

  const start = startOne - 1;
  const end = Math.min(endOne, length);
  if (end - start > maxBases) {
    throw new RegionError(
      `That region is ${(end - start).toLocaleString()} bases; the most you can extract ` +
        `in the browser is ${maxBases.toLocaleString()}. For bigger pieces, use one of ` +
        'the tools below.',
    );
  }
  return {
    name: seq.name,
    start,
    end,
    clamped: endOne > length,
    label: `${seq.name}:${start + 1}-${end}`,
  };
}

/** A FASTA record: `>label`, then the bases wrapped at `width`. */
export function regionFasta(label: string, bases: string, width = 80): string {
  const lines = [`>${label}`];
  for (let i = 0; i < bases.length; i += width) lines.push(bases.slice(i, i + width));
  return `${lines.join('\n')}\n`;
}

/** A file name for a region: `chr1_1000000-1000500.fa`. */
export function regionFileName(label: string): string {
  return `${label.replace(/[^\w.-]+/g, '_')}.fa`;
}

/** The default region for a sequence: its first 1,000 bases. */
export function defaultRegion(row: Pick<RegionSequence, 'name' | 'length'>): string {
  return `${row.name}:1-${Math.min(row.length, 1000)}`;
}
