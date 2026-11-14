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
  /** The name alone, for sentences. */
  name: string;
  /** Every alias, for type-to-filter. */
  aliases: string[];
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

function candidateName(genome: CandidateGenome): string {
  return preferredAlias(genome.aliases) ?? formatDigest(genome.digest);
}

export function candidateLabel(genome: CandidateGenome): string {
  const name = candidateName(genome);
  const detail = genome.assembly_accession || genome.assembly_source;
  return detail ? `${name} · ${detail}` : name;
}

/** An alias's leading letters and digits, lowercased: `hg38_primary` and `hg38-ddbj` give `hg38`. */
export function nameStem(alias: string): string | null {
  const match = /^[a-z]+\d+/i.exec(alias);
  return match ? match[0].toLowerCase() : null;
}

function stems(aliases: readonly string[]): Set<string> {
  return new Set(aliases.map(nameStem).filter((stem): stem is string => stem !== null));
}

/**
 * The other genomes of `self`'s species, ordered. First, by name, those whose
 * name shares a stem with any alias of `self` (for hg38: every `hg38_*`,
 * `GRCh38-*` and `hs38*`). Then the rest: same accession, same assembly (any version, GCA or
 * GCF), same source, then the others; by label within each group. `self`'s
 * aliases come from its own entry in `genomes` when it has none.
 */
export function rankCompareCandidates(
  self: GenomeFacts & { aliases?: readonly string[] },
  genomes: readonly CandidateGenome[],
): CompareCandidate[] {
  const selfAliases = self.aliases ?? genomes.find((g) => g.digest === self.digest)?.aliases ?? [];
  const selfStems = stems(selfAliases);
  return genomes
    .filter((g) => g.digest !== self.digest && sameSpecies(self, g))
    .map((g) => {
      // By the name shown, not every alias: some carry wrong ones (an hg19
      // with a GRCh38-* alias), and the group should look right at a glance.
      const name = candidateName(g);
      const stem = nameStem(name);
      return {
        candidate: { digest: g.digest, label: candidateLabel(g), name, aliases: g.aliases },
        // Within the family, by name alone: tiers would sink hg38_primary.
        rank: stem !== null && selfStems.has(stem) ? 0 : 1 + tier(self, g),
      };
    })
    .sort(
      (x, y) =>
        x.rank - y.rank ||
        x.candidate.name.localeCompare(y.candidate.name) ||
        x.candidate.label.localeCompare(y.candidate.label),
    )
    .map(({ candidate }) => candidate);
}
