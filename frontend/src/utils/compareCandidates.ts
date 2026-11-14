/** Which genomes a genome can sensibly be compared with, likeliest first. */

import { formatDigest } from './format';
import { preferredAlias } from './genomes';

export interface CandidateGenome {
  digest: string;
  aliases: string[];
  species_name: string | null;
  taxon_id: number | null;
  assembly_source: string | null;
  assembly_accession: string | null;
}

/** What ranking needs from the genome being compared. */
export type GenomeFacts = Omit<CandidateGenome, 'aliases'>;

export interface CompareCandidate {
  digest: string;
  /** The name and its accession or source, e.g. `hg38_primary · UCSC`. */
  label: string;
}

/** `GCA_000001405.15` -> `{ base: '000001405', version: '15' }`; GCA and GCF share a base. */
function parseAccession(
  accession: string | null,
): { base: string; version: string } | null {
  const match = accession?.trim().match(/^GC[AF]_(\d+)(?:\.(\d+))?$/i);
  return match ? { base: match[1], version: match[2] ?? '' } : null;
}

/** True if both genomes are the same species: by NCBI taxon id, else by species name. */
export function sameSpecies(a: GenomeFacts, b: GenomeFacts): boolean {
  if (a.taxon_id !== null && b.taxon_id !== null) return a.taxon_id === b.taxon_id;
  if (!a.species_name || !b.species_name) return false;
  return a.species_name.trim().toLowerCase() === b.species_name.trim().toLowerCase();
}

/** 0 = same accession, 1 = same assembly (accession base), 2 = same source, 3 = other. */
function tier(self: GenomeFacts, other: GenomeFacts): number {
  const a = parseAccession(self.assembly_accession);
  const b = parseAccession(other.assembly_accession);
  if (a && b && a.base === b.base) return a.version === b.version ? 0 : 1;
  const source = (g: GenomeFacts) => g.assembly_source?.trim().toLowerCase() || null;
  if (source(self) && source(self) === source(other)) return 2;
  return 3;
}

export function candidateLabel(genome: CandidateGenome): string {
  const name = preferredAlias(genome.aliases) ?? formatDigest(genome.digest);
  const detail = genome.assembly_accession || genome.assembly_source;
  return detail ? `${name} · ${detail}` : name;
}

/**
 * The other genomes of `self`'s species, ordered: same accession, same
 * assembly (any version, GCA or GCF), same source, then the rest; by label
 * within each group.
 */
export function rankCompareCandidates(
  self: GenomeFacts,
  genomes: readonly CandidateGenome[],
): CompareCandidate[] {
  return genomes
    .filter((g) => g.digest !== self.digest && sameSpecies(self, g))
    .map((g) => ({ digest: g.digest, label: candidateLabel(g), tier: tier(self, g) }))
    .sort((x, y) => x.tier - y.tier || x.label.localeCompare(y.label))
    .map(({ digest, label }) => ({ digest, label }));
}
