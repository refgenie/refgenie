/**
 * Haplotype from a genome's name. VGP-style ToLID names carry it as a token
 * (`fLatCha1.pri`, `bAegAlb1.hap1`, `aBomVar4.paternal.1`); there is no field.
 */

import type { GenomeResponse } from '../types/api';

export type Haplotype =
  'primary' | 'alternate' | 'hap1' | 'hap2' | 'maternal' | 'paternal';

export const HAPLOTYPES: readonly Haplotype[] = [
  'primary',
  'alternate',
  'hap1',
  'hap2',
  'maternal',
  'paternal',
] as const;

const TOKEN_TO_HAPLOTYPE: Record<string, Haplotype> = {
  pri: 'primary',
  primary: 'primary',
  alt: 'alternate',
  alternate: 'alternate',
  hap1: 'hap1',
  hap2: 'hap2',
  mat: 'maternal',
  maternal: 'maternal',
  pat: 'paternal',
  paternal: 'paternal',
};

/**
 * A delimiter or string edge on both sides, so a ToLID species code like
 * `mPatXyz1` never matches `pat`.
 */
const TOKEN =
  /(?:^|[._\s-])(pri|primary|alt|alternate|hap1|hap2|mat|maternal|pat|paternal)(?=$|[._\s-])/i;

function fromText(text: string | null | undefined): Haplotype | null {
  if (!text) return null;
  const match = TOKEN.exec(text);
  return match ? TOKEN_TO_HAPLOTYPE[match[1].toLowerCase()] : null;
}

/** Scan aliases, then synonyms, then the description; the first match wins. */
export function haplotypeOf(
  genome: Pick<GenomeResponse, 'aliases' | 'synonyms' | 'description'>,
): Haplotype | null {
  for (const text of [
    ...(genome.aliases ?? []),
    ...(genome.synonyms ?? []),
    genome.description,
  ]) {
    const found = fromText(text);
    if (found) return found;
  }
  return null;
}
