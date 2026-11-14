/** Pure helpers for the genome page's FASTA download: presets and the CLI command. */

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
 * bgzip bytes per FASTA base, as `bgzf.ts` writes it: measured 0.28 on hg38
 * and 0.30 on a genome with no N.
 */
export const BGZF_BYTES_PER_BASE = 0.3;

/** The likely size of a .fa.gz of `bases` bases. */
export function estimateFastaGzBytes(bases: number): number {
  return Math.round(bases * BGZF_BYTES_PER_BASE);
}
