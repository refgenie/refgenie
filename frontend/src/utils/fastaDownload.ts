/** Pure helpers for the genome page's FASTA download: presets and the CLI command. */
import { formatBytes } from './format';

/**
 * Names that are not primary chromosomes: unlocalized (_random), unplaced
 * (chrUn, Un_, GL/KI/JH contigs), alternate loci (_alt, _hap), patches
 * (_fix), decoys, EBV and HLA.
 */
const NON_PRIMARY =
  /(_random$|_alt$|_fix$|_hap\d+$|^chrUn|^Un[_.]|_Un_|decoy|^(chr)?EBV$|^HLA[-_]|^(GL|KI|JH|KN|KV|KZ|ML|MU)\d+)/i;

/**
 * The primary chromosomes, by name. `null` when the pattern does not split
 * the genome (it keeps all names or none), so offering it would not help.
 */
export function primaryChromosomes(names: string[]): string[] | null {
  const primary = names.filter((name) => !NON_PRIMARY.test(name));
  if (primary.length === 0 || primary.length === names.length) return null;
  return primary;
}

/** Quote a word for a POSIX shell, only when it needs it. */
export function shellQuote(word: string): string {
  if (/^[A-Za-z0-9._:/@%+=,#-]+$/.test(word)) return word;
  return `'${word.replace(/'/g, `'\\''`)}'`;
}

/**
 * The `refget store export` command (refget 0.13) that writes the same
 * sequences the browser does. `names` limits it to those sequences; omit it
 * for all of them. Unlike the browser file, whose headers carry names only,
 * its headers also keep each sequence's description.
 */
export function refgetExportCommand({
  digest,
  storeUrl,
  output,
  names,
}: {
  digest: string;
  storeUrl: string;
  output: string;
  names?: string[];
}): string {
  const words = [
    'refget',
    'store',
    'export',
    digest,
    '--remote',
    storeUrl,
    ...(names ?? []).flatMap((name) => ['-n', name]),
    '-o',
    output,
  ];
  return words.map(shellQuote).join(' ');
}

/**
 * bgzip bytes per base of FASTA (80-column lines, uppercase), as written by
 * `bgzf.ts`. Measured on 48 evenly spread 256 kb windows of two genomes:
 * hg38 gave 0.281 with 6% N (about 0.30 per non-N base; N runs cost almost
 * nothing), and a VGP bat assembly with no N gave 0.303. The sequence list
 * does not say how much of a sequence is N, so every base counts as 0.30; for
 * hg38 that overshoots by about 5%.
 */
export const BGZF_BYTES_PER_BASE = 0.3;

/** The likely size of a .fa.gz of `bases` bases. */
export function estimateFastaGzBytes(bases: number): number {
  return Math.round(bases * BGZF_BYTES_PER_BASE);
}

/**
 * A rough size with two significant figures: "0.9 GB", "640 MB", "12 KB".
 * Anything from half a unit up is shown in that unit ("0.6 GB", not "610 MB").
 */
export function approxBytes(bytes: number): string {
  const units = ['B', 'KB', 'MB', 'GB', 'TB'];
  let value = bytes;
  let unit = 0;
  while (value >= 512 && unit < units.length - 1) {
    value /= 1024;
    unit += 1;
  }
  if (unit === 0) return `${Math.round(value)} B`;
  const shown = value < 10 ? value.toFixed(1) : String(Math.round(value));
  return `${shown} ${units[unit]}`;
}

/** "Saving… 34% · 312.0 MB of about 920.0 MB": the total is projected from progress so far. */
export function savingText(fraction: number, written: number): string {
  const percent = Math.floor(fraction * 100);
  const total = fraction >= 0.01 ? ` of about ${formatBytes(written / fraction)}` : '';
  return `Saving… ${percent}% · ${formatBytes(written)}${total}`;
}
