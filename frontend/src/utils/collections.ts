/** Pure helpers for collection pages: headline counts, per-genome links, filtering. */

import { haplotypeOf } from './haplotype';
import type { Haplotype } from './haplotype';
import { normalizeSpecies } from './species';
import type { GenomeResponse } from '../types/api';

export interface CollectionStats {
  genomes: number;
  species: number;
  withAssets: number;
  assetGroups: number;
  haplotypes: Record<Haplotype | 'unspecified', number>;
}

export function collectionStats(genomes: readonly GenomeResponse[]): CollectionStats {
  const species = new Set<string>();
  const haplotypes: CollectionStats['haplotypes'] = {
    primary: 0,
    alternate: 0,
    hap1: 0,
    hap2: 0,
    maternal: 0,
    paternal: 0,
    unspecified: 0,
  };
  let withAssets = 0;
  let assetGroups = 0;
  for (const genome of genomes) {
    if (genome.species_name) species.add(normalizeSpecies(genome.species_name));
    if (genome.asset_count > 0) withAssets += 1;
    assetGroups += genome.asset_count;
    haplotypes[haplotypeOf(genome) ?? 'unspecified'] += 1;
  }
  return {
    genomes: genomes.length,
    species: species.size,
    withAssets,
    assetGroups,
    haplotypes,
  };
}

type LinkGenome = Pick<
  GenomeResponse,
  'assembly_accession' | 'species_name' | 'digest'
>;

function placeholderValue(name: string, genome: LinkGenome): string | null | undefined {
  switch (name) {
    case 'assembly_accession':
      return genome.assembly_accession;
    case 'species_name':
      return genome.species_name;
    case 'species_underscored':
      return genome.species_name?.trim().replace(/\s+/g, '_');
    case 'digest':
      return genome.digest;
    default:
      return undefined;
  }
}

/**
 * Fill a per-genome link template. Returns `null` if any placeholder it uses
 * is empty or unknown, so no half-filled link is ever rendered.
 */
export function expandGenomeLink(template: string, genome: LinkGenome): string | null {
  let ok = true;
  const url = template.replace(/\{([a-z_]+)\}/g, (_, name: string) => {
    const value = placeholderValue(name, genome);
    if (!value) {
      ok = false;
      return '';
    }
    return encodeURIComponent(value);
  });
  return ok ? url : null;
}

/** True when the template fills in `placeholder`. */
export function usesPlaceholder(template: string, placeholder: string): boolean {
  return template.includes(`{${placeholder}}`);
}

export interface CollectionFilter {
  q?: string;
  haplotype?: Haplotype | null;
  onlyWithAssets?: boolean;
}

/** Client-side filter: a case-insensitive substring over the genome's names and accession. */
export function filterCollectionGenomes(
  genomes: readonly GenomeResponse[],
  { q, haplotype, onlyWithAssets }: CollectionFilter,
): GenomeResponse[] {
  const needle = q?.trim().toLowerCase() ?? '';
  return genomes.filter((genome) => {
    if (onlyWithAssets && genome.asset_count <= 0) return false;
    if (haplotype && haplotypeOf(genome) !== haplotype) return false;
    if (!needle) return true;
    const haystack = [
      ...genome.aliases,
      genome.species_name,
      genome.common_name,
      genome.assembly_accession,
      genome.description,
    ];
    return haystack.some((text) => text?.toLowerCase().includes(needle));
  });
}
