/** Digest -> genome. Built once per session by `useGenomeIndex()`. */

import { formatDigest } from './format';
import type { GenomeResponse } from '../types/api';
import type { SearchOperator } from '../types/pagination';

export type GenomeIndex = Map<string, GenomeResponse>;

/**
 * Module-level so it can be passed to the cache's `select` by reference. The
 * memo is keyed on this function's identity: an inline arrow would be a new
 * function every render and the Map would be rebuilt every render with it.
 */
export function buildGenomeIndex(genomes: readonly GenomeResponse[]): GenomeIndex {
  return new Map(genomes.map((genome) => [genome.digest, genome]));
}

/**
 * The alias to lead with.
 *
 * `aliases[0]` is insertion order, not canonicality: the public server's human
 * genome carries
 * `['hg38-refgenie', 'hg38-noALT-noHLA-noDecoy-broad', 'hg38', …]`, and a row
 * matched by a search for `hg38` should echo the term that found it, not label
 * itself `hg38-refgenie`.
 *
 * Prefer an exact match on what the user typed, then a substring match, then
 * the shortest alias — short aliases are the canonical assembly names
 * (`hg38`, `mm10`) and long ones are qualified variants of them.
 */
export function preferredAlias(
  aliases: readonly string[] | null | undefined,
  query?: string,
): string | undefined {
  if (!aliases || aliases.length === 0) return undefined;
  const needle = query?.trim().toLowerCase();
  if (needle) {
    const exact = aliases.find((alias) => alias.toLowerCase() === needle);
    if (exact) return exact;
    const partial = aliases.find((alias) => alias.toLowerCase().includes(needle));
    if (partial) return partial;
  }
  return [...aliases].sort((a, b) => a.length - b.length || a.localeCompare(b))[0];
}

/**
 * The label a human recognizes, falling back to a truncated digest so a row is
 * never blank while the index loads or when a genome carries no alias.
 */
export function genomeLabel(
  index: GenomeIndex | undefined,
  digest: string | null | undefined,
  query?: string,
): string {
  if (!digest) return 'NA';
  return preferredAlias(index?.get(digest)?.aliases, query) ?? formatDigest(digest);
}

/**
 * The searchable text of one genome, keyed by the field names
 * `GENOME_SEARCH_FIELDS` declares. Mirrors `catalog.py::list_genomes`.
 */
const GENOME_FIELD_VALUES: Record<string, (genome: GenomeResponse) => string[]> = {
  digest: (genome) => [genome.digest],
  description: (genome) => (genome.description ? [genome.description] : []),
  species_name: (genome) => (genome.species_name ? [genome.species_name] : []),
  common_name: (genome) => (genome.common_name ? [genome.common_name] : []),
  assembly_source: (genome) => (genome.assembly_source ? [genome.assembly_source] : []),
  assembly_accession: (genome) =>
    genome.assembly_accession ? [genome.assembly_accession] : [],
  aliases: (genome) => [...genome.aliases],
};

/** `eq` is a case-sensitive `==`; the rest are case-insensitive, like SQL `ilike`. */
function matchesTerm(value: string, needle: string, operator: SearchOperator): boolean {
  if (operator === 'eq') return value === needle;
  const haystack = value.toLowerCase();
  const term = needle.toLowerCase();
  if (operator === 'starts_with') return haystack.startsWith(term);
  if (operator === 'ends_with') return haystack.endsWith(term);
  return haystack.includes(term);
}

/**
 * The client-side twin of the server's genome search.
 *
 * `GenomesPage` filters in memory when its "only genomes with assets" filter is
 * on, because no endpoint has a `has_assets` parameter. This keeps the search
 * box, its field selector and its operator selector meaning the same thing in
 * both modes; without it they would silently do nothing while filtering.
 */
export function genomeMatches(
  genome: GenomeResponse,
  query: string,
  fields: readonly string[],
  operator: SearchOperator,
): boolean {
  const needle = query.trim();
  if (!needle) return true;
  const selected = fields.length > 0 ? fields : Object.keys(GENOME_FIELD_VALUES);
  return selected
    .flatMap((field) => GENOME_FIELD_VALUES[field]?.(genome) ?? [])
    .some((value) => matchesTerm(value, needle, operator));
}
