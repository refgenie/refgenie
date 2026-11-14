/** Files used with a genome must name its sequences the same way, so this is worth a badge. */
export const NAMING_TITLE =
  'Files you use with this genome (BED, GTF, VCF) must use the same sequence names.';

/**
 * The naming style of a genome's longest sequence: `UCSC-style names (chr1)`
 * for a `chr` prefix, `Ensembl-style names (1)` for a short plain name, else null.
 */
export function namingStyle(
  sequences: ReadonlyArray<{ name: string; length: number }>,
): string | null {
  if (sequences.length === 0) return null;
  const { name } = sequences.reduce((a, b) => (b.length > a.length ? b : a));
  if (/^chr/i.test(name)) return `UCSC-style names (${name})`;
  return /^[A-Z0-9]{1,4}$/i.test(name) ? `Ensembl-style names (${name})` : null;
}
