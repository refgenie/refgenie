/**
 * Which naming convention a genome's sequences follow: UCSC's `chr1` and
 * `chrM`, Ensembl's `1` and `MT`, or GenBank/RefSeq accessions. Files used
 * with a genome must name sequences the same way, so this is worth a badge.
 */

export type NamingStyle = 'ucsc' | 'ensembl' | 'accession';

export interface Naming {
  style: NamingStyle;
  /** Names from this genome that show the style, e.g. `['chr1', 'chrM']`. */
  examples: string[];
}

const MITO = /^(chr)?(M|MT|Mito)$/i;
/** `NC_000001.11`, `CM000663.2`, `KI270728.1`. */
const ACCESSION = /^[A-Z]{1,2}_?\d{5,}(\.\d+)?$/;
/** `1`, `X`, `MT`, `2L`, `III`: short, no prefix. */
const SHORT = /^[A-Za-z0-9]{1,4}$/;

function styleOf(name: string): NamingStyle | null {
  if (/^chr/i.test(name)) return 'ucsc';
  if (ACCESSION.test(name)) return 'accession';
  if (SHORT.test(name)) return 'ensembl';
  return null;
}

/**
 * The naming style that covers most of a genome's bases (so its chromosomes
 * decide, not its many small scaffolds), or null when no style covers more
 * than half.
 */
export function namingStyle(
  sequences: ReadonlyArray<{ name: string; length: number }>,
): Naming | null {
  const bases = new Map<NamingStyle, number>();
  const longest = new Map<NamingStyle, { name: string; length: number }>();
  let total = 0;
  for (const s of sequences) {
    total += s.length;
    const style = styleOf(s.name);
    if (!style) continue;
    bases.set(style, (bases.get(style) ?? 0) + s.length);
    if (!MITO.test(s.name) && s.length > (longest.get(style)?.length ?? -1)) {
      longest.set(style, s);
    }
  }
  const [style, covered] = [...bases.entries()].sort((a, b) => b[1] - a[1])[0] ?? [];
  if (!style || covered === undefined || covered * 2 <= total) return null;
  const mito = sequences.find((s) => MITO.test(s.name) && styleOf(s.name) === style);
  const examples = [longest.get(style)?.name, mito?.name].filter((n): n is string => !!n);
  return { style, examples };
}

const LABELS: Record<NamingStyle, string> = {
  ucsc: 'UCSC-style names',
  ensembl: 'Ensembl-style names',
  accession: 'Accession names',
};

const TITLES: Record<NamingStyle, string> = {
  ucsc: 'Sequence names start with "chr", as at UCSC. Files you use with this genome (BED, GTF, VCF) must use the same names: chr1, not 1.',
  ensembl:
    'Sequence names have no "chr" prefix, as at Ensembl. Files you use with this genome (BED, GTF, VCF) must use the same names: 1, not chr1.',
  accession:
    'Sequences are named by their GenBank or RefSeq accession. Files you use with this genome must use the same names.',
};

/** `UCSC-style names (chr1, chrM)`. */
export function namingLabel(naming: Naming): string {
  const label = LABELS[naming.style];
  return naming.examples.length ? `${label} (${naming.examples.join(', ')})` : label;
}

export function namingTitle(naming: Naming): string {
  return TITLES[naming.style];
}
