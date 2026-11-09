/**
 * Read a BED file of regions to extract from a genome in the browser.
 *
 * BED coordinates are 0-based and half-open: `chr1 999999 1000500` is the 501
 * bases from position 1,000,000 to 1,000,500 (1-based). Columns may be split
 * by tabs or spaces; `track`, `browser` and `#` lines are skipped; a 4th
 * column is the region's name. Strand (column 6) is ignored: bases are
 * always read from the + strand.
 *
 * Each region becomes one FASTA record, labeled as `bedtools getfasta -name`
 * does: `name::chr:start-end` with a name column, else `chr:start-end`, with
 * the BED's own coordinates and sequence name.
 */
import {
  exampleRegion,
  findSequence,
  suggestSequences,
  type RegionSequence,
} from './parseRegion';

/** The most regions one in-browser BED extraction may read. */
export const MAX_BED_REGIONS = 5_000;
/** The most bases, summed over all regions, one in-browser BED extraction may read. */
export const MAX_BED_BASES = 5_000_000;
/** The largest BED file (in bytes, after unzipping) the browser will read. */
export const MAX_BED_FILE_BYTES = 20_000_000;

export interface BedSequence extends RegionSequence {
  alphabet?: string | null;
}

export interface BedRegion {
  /** 1-based line number in the BED text. */
  line: number;
  /** The sequence name as written in the BED (a name or an alias). */
  chrom: string;
  /** The genome's own name for that sequence. */
  name: string;
  /** 0-based, inclusive. */
  start: number;
  /** 0-based, exclusive. */
  end: number;
  /** The BED name column, if any. */
  label: string | null;
  /** The FASTA header, without the `>`. */
  header: string;
}

export interface BedProblem {
  line: number;
  text: string;
  message: string;
}

export interface ParsedBed {
  regions: BedRegion[];
  problems: BedProblem[];
  /** Total bases over all good regions. */
  bases: number;
  /** Header, comment and blank lines skipped. */
  skipped: number;
}

const SKIP = /^(#|track(\s|$)|browser(\s|$))/;
const INTEGER = /^\d+$/;
const fmt = (n: number) => n.toLocaleString('en-US');

/** The FASTA header for a BED region, as `bedtools getfasta -name` writes it. */
export function bedHeader(
  chrom: string,
  start: number,
  end: number,
  label?: string | null,
): string {
  const coords = `${chrom}:${start}-${end}`;
  return label ? `${label}::${coords}` : coords;
}

/**
 * Parse BED text against a genome's sequences. Lines that cannot be read, name
 * an unknown sequence, or fall outside it are returned as problems; the rest
 * are regions ready to read.
 */
export function parseBed(
  text: string,
  sequences: readonly BedSequence[],
  isReadable: (seq: BedSequence) => boolean = () => true,
): ParsedBed {
  const regions: BedRegion[] = [];
  const problems: BedProblem[] = [];
  let bases = 0;
  let skipped = 0;
  const lines = text.split('\n');
  // Name lookups repeat for every line on the same sequence.
  const found = new Map<string, BedSequence | null>();
  const lookup = (chrom: string) => {
    if (!found.has(chrom)) found.set(chrom, findSequence(chrom, sequences));
    return found.get(chrom) ?? null;
  };

  for (let i = 0; i < lines.length; i += 1) {
    const raw = lines[i].replace(/\r$/, '');
    const trimmed = raw.trim();
    if (trimmed === '' || SKIP.test(trimmed)) {
      // The empty "line" after a final newline is not a skipped line.
      if (trimmed !== '' || i < lines.length - 1) skipped += 1;
      continue;
    }
    const line = i + 1;
    const problem = (message: string) => problems.push({ line, text: raw, message });
    // Tabs split columns when present, so a name may contain spaces.
    const cols = (raw.includes('\t') ? raw.split('\t') : trimmed.split(/\s+/)).map(
      (c) => c.trim(),
    );
    const [chrom, startText, endText, labelText] = cols;
    if (!chrom || !INTEGER.test(startText ?? '') || !INTEGER.test(endText ?? '')) {
      problem('Not a BED line: it needs a sequence name, a start and an end.');
      continue;
    }
    const start = Number(startText);
    const end = Number(endText);
    const seq = lookup(chrom);
    if (!seq) {
      const close = suggestSequences(chrom, sequences);
      problem(
        `No sequence named "${chrom}".` +
          (close.length
            ? ` Did you mean ${close.map((n) => `"${n}"`).join(' or ')}?`
            : ''),
      );
      continue;
    }
    if (start >= end) {
      problem(`The start (${fmt(start)}) must be less than the end (${fmt(end)}).`);
      continue;
    }
    if (end > seq.length) {
      problem(
        `The end (${fmt(end)}) is past the end of ${seq.name} (${fmt(seq.length)} bases).`,
      );
      continue;
    }
    if (!isReadable(seq)) {
      problem(
        `${seq.name} uses the ${seq.alphabet ?? 'unknown'} alphabet, which cannot be read here.`,
      );
      continue;
    }
    const label = labelText && labelText !== '.' ? labelText : null;
    regions.push({
      line,
      chrom,
      name: seq.name,
      start,
      end,
      label,
      header: bedHeader(chrom, start, end, label),
    });
    bases += end - start;
  }
  return { regions, problems, bases, skipped };
}

/** Why a parsed BED is too big to run in the browser, or null if it fits. */
export function bedOverCap(
  bed: Pick<ParsedBed, 'regions' | 'bases'>,
  maxRegions = MAX_BED_REGIONS,
  maxBases = MAX_BED_BASES,
): string | null {
  const parts: string[] = [];
  if (bed.regions.length > maxRegions) {
    parts.push(`${fmt(bed.regions.length)} regions (the limit is ${fmt(maxRegions)})`);
  }
  if (bed.bases > maxBases) {
    parts.push(`${fmt(bed.bases)} bases (the limit is ${fmt(maxBases)})`);
  }
  return parts.length ? `This file has ${parts.join(' and ')}.` : null;
}

/** One FASTA record per region, each sequence on one line, as bedtools getfasta writes. */
export function bedFastaRecord(header: string, bases: string): string {
  return `>${header}\n${bases}\n`;
}

/**
 * A small BED file for this genome: three named regions on its longest
 * sequences (the same spot as the single-region example on the longest), and
 * one unnamed region at the start of the shortest.
 */
export function exampleBed(sequences: readonly RegionSequence[]): string | null {
  const ranked = [...sequences]
    .filter((s) => s.length > 0)
    .sort((a, b) => b.length - a.length);
  if (!ranked.length) return null;
  const lines: string[] = [];
  const one = exampleRegion(ranked);
  if (one) {
    // "name:1,000,000-1,000,500" (1-based) is BED 999999 1000500.
    const match = /^(.*):([\d,]+)-([\d,]+)$/.exec(one);
    if (match) {
      const start = Number(match[2].replace(/,/g, '')) - 1;
      const end = Number(match[3].replace(/,/g, ''));
      lines.push([match[1], start, end, 'region_1'].join('\t'));
    }
  }
  ranked.slice(1, 3).forEach((seq, i) => {
    const start = seq.length > 2_000_000 ? 2_000_000 : Math.floor(seq.length / 2);
    const end = Math.min(seq.length, start + 300);
    lines.push([seq.name, start, end, `region_${i + 2}`].join('\t'));
  });
  const shortest = ranked[ranked.length - 1];
  if (ranked.length > 1)
    lines.push([shortest.name, 0, Math.min(shortest.length, 60)].join('\t'));
  return `${lines.join('\n')}\n`;
}

/** A FASTA file name for a BED file: `peaks.bed.gz` becomes `peaks.fa`. */
export function bedFastaFileName(bedName: string | null | undefined): string {
  const stem = (bedName ?? '').replace(/\.gz$/i, '').replace(/\.(bed|txt|tsv)$/i, '');
  return `${stem.replace(/[^\w.-]+/g, '_') || 'regions'}.fa`;
}

/** True if bytes start with the gzip magic number. */
export function isGzip(bytes: Uint8Array): boolean {
  return bytes.length >= 2 && bytes[0] === 0x1f && bytes[1] === 0x8b;
}

/** A BED file's text, unzipping a .gz file (by its first bytes, not its name). */
export async function readBedFile(file: Blob): Promise<string> {
  if (file.size > MAX_BED_FILE_BYTES * 4) {
    throw new Error(`That file is ${fmt(file.size)} bytes, too big to read here.`);
  }
  const bytes = await blobBytes(file);
  if (!isGzip(bytes)) {
    if (bytes.length > MAX_BED_FILE_BYTES) {
      throw new Error(`That file is ${fmt(bytes.length)} bytes, too big to read here.`);
    }
    return new TextDecoder().decode(bytes);
  }
  if (typeof DecompressionStream === 'undefined') {
    throw new Error('This browser cannot unzip .gz files; unzip it first.');
  }
  const zipped = new Response(bytes as BlobPart).body as ReadableStream<BufferSource>;
  const text = await new Response(
    zipped.pipeThrough(new DecompressionStream('gzip')),
  ).text();
  if (text.length > MAX_BED_FILE_BYTES) {
    throw new Error('That file is too big to read here once unzipped.');
  }
  return text;
}

/** A Blob's bytes, by FileReader where Blob.arrayBuffer is missing (jsdom). */
function blobBytes(blob: Blob): Promise<Uint8Array> {
  if (typeof blob.arrayBuffer === 'function') {
    return blob.arrayBuffer().then((buffer) => new Uint8Array(buffer));
  }
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(new Uint8Array(reader.result as ArrayBuffer));
    reader.onerror = () => reject(reader.error);
    reader.readAsArrayBuffer(blob);
  });
}
