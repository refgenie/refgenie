/**
 * Species rows grouped in the browser from the cached genome list. We skip
 * `/v4/species/summary` because it lacks common names, drops species with no
 * assets, and is not served by a local dash.
 */

import type { GenomeResponse } from '../types/api';

export interface SpeciesRow {
  /** Grouping key: lowercased, whitespace-collapsed. '' = no species recorded. */
  key: string;
  /** The spelling to render. '' rows render as "Unspecified". */
  speciesName: string;
  /** Every distinct non-null `common_name` seen for this species, sorted. */
  commonNames: string[];
  /** Other names from the genomes' `species` objects; for search only. */
  otherNames: string[];
  /** Usually one; more than one is a data problem and is shown. */
  taxonIds: number[];
  genomes: number;
  assets: number;
}

export type SpeciesIndex = Map<string, SpeciesRow>;

/**
 * The grouping key, case-folded because `taxa.json` and the API spell some
 * names differently ("Homo Sapiens" vs "Homo sapiens").
 */
export function normalizeSpecies(name: string | null | undefined): string {
  return (name ?? '').trim().replace(/\s+/g, ' ').toLowerCase();
}

/** Module-level so it can be handed to the cache's `select` by reference. */
export function buildSpeciesIndex(genomes: readonly GenomeResponse[]): SpeciesIndex {
  // Count spellings so the row shows the most common one.
  const spellings = new Map<string, Map<string, number>>();
  const commons = new Map<string, Set<string>>();
  const others = new Map<string, Set<string>>();
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
        otherNames: [],
        taxonIds: [],
        genomes: 1,
        assets: genome.asset_count,
      });
      spellings.set(key, new Map());
      commons.set(key, new Set());
      others.set(key, new Set());
      taxa.set(key, new Set());
    }

    if (genome.species_name) {
      const counts = spellings.get(key)!;
      counts.set(genome.species_name, (counts.get(genome.species_name) ?? 0) + 1);
    }
    if (genome.common_name) commons.get(key)!.add(genome.common_name);
    if (genome.species) {
      const { scientific_name, common_name, synonyms } = genome.species;
      for (const name of [scientific_name, common_name, ...synonyms]) {
        if (name) others.get(key)!.add(name);
      }
    }
    if (genome.taxon_id !== null) taxa.get(key)!.add(genome.taxon_id);
  }

  for (const [key, row] of index) {
    const counts = [...spellings.get(key)!.entries()];
    // Majority spelling; alphabetical break so the render is deterministic.
    counts.sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
    row.speciesName = counts[0]?.[0] ?? '';
    row.commonNames = [...commons.get(key)!].sort((a, b) => a.localeCompare(b));
    row.otherNames = [...others.get(key)!].sort((a, b) => a.localeCompare(b));
    row.taxonIds = [...taxa.get(key)!].sort((a, b) => a - b);
  }

  return index;
}

/**
 * Case-insensitive substring match on the scientific name, common names, other
 * names, and taxon ID (users search for `9606`).
 */
export function speciesMatches(row: SpeciesRow, query: string): boolean {
  const needle = query.trim().toLowerCase();
  if (!needle) return true;
  if (row.key.includes(needle)) return true;
  if (row.commonNames.some((name) => name.toLowerCase().includes(needle))) return true;
  if (row.otherNames.some((name) => name.toLowerCase().includes(needle))) return true;
  return row.taxonIds.some((id) => String(id).includes(needle));
}

/** Unspecified last; then species with the most built assets. */
export function compareSpecies(a: SpeciesRow, b: SpeciesRow): number {
  if ((a.key === '') !== (b.key === '')) return a.key === '' ? 1 : -1;
  return (
    b.assets - a.assets ||
    b.genomes - a.genomes ||
    a.speciesName.localeCompare(b.speciesName)
  );
}

/**
 * Link to the genome list filtered to this species, or null for the unspecified
 * row. Uses `op=contains` because server `eq` is case-sensitive and rows are
 * grouped case-insensitively.
 */
export function speciesGenomesHref(row: SpeciesRow): string | null {
  if (row.key === '') return null;
  const params = new URLSearchParams({
    q: row.speciesName,
    fields: 'species',
    op: 'contains',
  });
  return `/genomes?${params.toString()}`;
}
