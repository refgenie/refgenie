/**
 * Species rows derived from the genome list.
 *
 * Nothing here is a wire type: `/v4/species/summary` is not called by this app
 * (it carries no common name, omits every species whose genomes have no assets,
 * and is not mounted on a local dash). Species are grouped in the browser from
 * the already-cached `useGenomeIndex()` payload instead.
 */

import type { GenomeResponse } from '../types/api';

export interface SpeciesRow {
  /** Grouping key: lowercased, whitespace-collapsed. '' = no species recorded. */
  key: string;
  /** The spelling to render. '' rows render as "Unspecified". */
  speciesName: string;
  /** Every distinct non-null `common_name` seen for this species, sorted. */
  commonNames: string[];
  /** Usually one. More than one is a data problem, so it is shown, not hidden. */
  taxonIds: number[];
  genomes: number;
  assets: number;
}

export type SpeciesIndex = Map<string, SpeciesRow>;

/**
 * The grouping key.
 *
 * Case-folded on purpose: `taxa.json` says "Homo Sapiens" and the API says
 * "Homo sapiens", and without folding the tree reports no genomes for humans
 * because of exactly that mismatch.
 */
export function normalizeSpecies(name: string | null | undefined): string {
  return (name ?? '').trim().replace(/\s+/g, ' ').toLowerCase();
}

/** Module-level so it can be handed to the cache's `select` by reference. */
export function buildSpeciesIndex(genomes: readonly GenomeResponse[]): SpeciesIndex {
  // Spelling counts, so the row renders the majority spelling rather than
  // whichever row happened to be first out of the database.
  const spellings = new Map<string, Map<string, number>>();
  const commons = new Map<string, Set<string>>();
  const taxa = new Map<string, Set<number>>();
  const index: SpeciesIndex = new Map();

  for (const genome of genomes) {
    const key = normalizeSpecies(genome.species_name);
    const row = index.get(key);
    if (row) {
      row.genomes += 1;
      row.assets += genome.asset_count;
    } else {
      index.set(key, {
        key,
        speciesName: '',
        commonNames: [],
        taxonIds: [],
        genomes: 1,
        assets: genome.asset_count,
      });
      spellings.set(key, new Map());
      commons.set(key, new Set());
      taxa.set(key, new Set());
    }

    if (genome.species_name) {
      const counts = spellings.get(key)!;
      counts.set(genome.species_name, (counts.get(genome.species_name) ?? 0) + 1);
    }
    if (genome.common_name) commons.get(key)!.add(genome.common_name);
    if (genome.taxon_id !== null) taxa.get(key)!.add(genome.taxon_id);
  }

  for (const [key, row] of index) {
    const counts = [...spellings.get(key)!.entries()];
    // Majority spelling; alphabetical break so the render is deterministic.
    counts.sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
    row.speciesName = counts[0]?.[0] ?? '';
    row.commonNames = [...commons.get(key)!].sort((a, b) => a.localeCompare(b));
    row.taxonIds = [...taxa.get(key)!].sort((a, b) => a - b);
  }

  return index;
}

/**
 * One box, both name kinds.
 *
 * Case-insensitive substring over the scientific name, every common name, and
 * the taxon ID — a bare `9606` is a search a genomics user will try.
 */
export function speciesMatches(row: SpeciesRow, query: string): boolean {
  const needle = query.trim().toLowerCase();
  if (!needle) return true;
  if (row.key.includes(needle)) return true;
  if (row.commonNames.some((name) => name.toLowerCase().includes(needle))) return true;
  return row.taxonIds.some((id) => String(id).includes(needle));
}

/** Unspecified last; then the species you can actually get files for. */
export function compareSpecies(a: SpeciesRow, b: SpeciesRow): number {
  if ((a.key === '') !== (b.key === '')) return a.key === '' ? 1 : -1;
  return (
    b.assets - a.assets ||
    b.genomes - a.genomes ||
    a.speciesName.localeCompare(b.speciesName)
  );
}

/**
 * The genome list filtered to this species, or null when no query can express
 * the row (`species_name IS NULL` is not a search term).
 *
 * `op=contains`, not `op=eq`: rows are grouped case-insensitively, and `eq` is
 * a case-sensitive `==` on both the server and in `genomeMatches`, so `eq`
 * would drop every genome spelled differently from the majority. The cost is
 * that a species name also matches its own subspecies, which is the harmless
 * direction. `assets=all` is added only when the species has nothing built, so
 * a click never lands on GenomesPage's "untick the filter" empty state.
 */
export function speciesGenomesHref(row: SpeciesRow): string | null {
  if (row.key === '') return null;
  const params = new URLSearchParams({
    q: row.speciesName,
    fields: 'species_name',
    op: 'contains',
  });
  if (row.assets === 0) params.set('assets', 'all');
  return `/genomes?${params.toString()}`;
}
