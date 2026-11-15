/** Digest -> genome. Built once per session by `useGenomeIndex()`. */

import { formatDigest } from './format';
import type { GenomeResponse } from '../types/api';

export type GenomeIndex = Map<string, GenomeResponse>;

/** Module-level so the cache's `select` memo sees a stable function. */
export function buildGenomeIndex(genomes: readonly GenomeResponse[]): GenomeIndex {
  return new Map(genomes.map((genome) => [genome.digest, genome]));
}

/**
 * The alias to show: an exact match on the query, then a substring match, then
 * the shortest alias that is not one of the genome's FHR `synonyms` (FHR lists
 * every name but the registry's own as a synonym, so `hg38_initial` wins over
 * `hg38-ddbj`), then the shortest alias.
 */
export function preferredAlias(
  aliases: readonly string[] | null | undefined,
  query?: string,
  synonyms?: readonly string[] | null,
): string | undefined {
  if (!aliases || aliases.length === 0) return undefined;
  const needle = query?.trim().toLowerCase();
  if (needle) {
    const exact = aliases.find((alias) => alias.toLowerCase() === needle);
    if (exact) return exact;
    const partial = aliases.find((alias) => alias.toLowerCase().includes(needle));
    if (partial) return partial;
  }
  const own = aliases.filter((alias) => !synonyms?.includes(alias));
  return [...(own.length ? own : aliases)].sort(
    (a, b) => a.length - b.length || a.localeCompare(b),
  )[0];
}

/** The genome's preferred alias, or a short digest when it has none. */
export function genomeLabel(
  index: GenomeIndex | undefined,
  digest: string | null | undefined,
  query?: string,
): string {
  if (!digest) return 'NA';
  return preferredAlias(index?.get(digest)?.aliases, query) ?? formatDigest(digest);
}
