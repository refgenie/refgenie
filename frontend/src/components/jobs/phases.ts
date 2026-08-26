/**
 * Phase vocabulary.
 *
 * The strings come from the JobManager and mirror the real phases of the
 * underlying code, so the card can say "step 4 of 8" instead of showing an
 * anonymous spinner.
 *
 * Three of them — `verify`, `digest`, `stage` — are full sha256 passes over
 * multi-gigabyte data that emit NOTHING while they run. They look exactly like
 * a hang, so each carries explicit copy saying it is normal and slow.
 */

import type { JobKind } from '../../services/contracts';

export interface PhaseInfo {
  label: string;
  /** Copy for the phases that produce no output for minutes at a time. */
  slowNote?: string;
}

const PULL_PHASES: Record<string, PhaseInfo> = {
  resolve: { label: 'Resolving genome' },
  query: { label: 'Querying server' },
  stage_lookup: { label: 'Looking up archive' },
  download: { label: 'Downloading' },
  verify: {
    label: 'Verifying checksum',
    slowNote: 'This can take several minutes for large assets, with no output.',
  },
  unpack: { label: 'Unpacking' },
  register: { label: 'Registering asset' },
  symlink: { label: 'Linking aliases' },
};

const BUILD_PHASES: Record<string, PhaseInfo> = {
  resolve_recipe: { label: 'Resolving recipe' },
  seek_keys: { label: 'Resolving seek keys' },
  validate: { label: 'Validating inputs' },
  run: { label: 'Running build' },
  digest: {
    label: 'Digesting output',
    slowNote: 'Hashing every file in the asset. Silent, and slow on large assets.',
  },
  register: { label: 'Registering asset' },
  stage: {
    label: 'Staging',
    slowNote: 'Creating a tar archive and hashing it. Silent, and slow on large assets.',
  },
};

const GENOME_INIT_PHASES: Record<string, PhaseInfo> = {
  resolve: { label: 'Reading FASTA' },
  digest: {
    label: 'Digesting sequences',
    slowNote: 'Hashing every sequence. Silent, and slow on large genomes.',
  },
  register: { label: 'Registering genome' },
  build: { label: 'Building fasta asset' },
};

export const PHASE_ORDER: Record<JobKind, readonly string[]> = {
  pull: ['resolve', 'query', 'stage_lookup', 'download', 'verify', 'unpack', 'register', 'symlink'],
  build: ['resolve_recipe', 'seek_keys', 'validate', 'run', 'digest', 'register', 'stage'],
  genome_init: ['resolve', 'digest', 'register', 'build'],
};

const BY_KIND: Record<JobKind, Record<string, PhaseInfo>> = {
  pull: PULL_PHASES,
  build: BUILD_PHASES,
  genome_init: GENOME_INIT_PHASES,
};

function humanize(phase: string): string {
  return phase.replace(/_/g, ' ').replace(/^./, (c) => c.toUpperCase());
}

export function phaseInfo(kind: JobKind, phase: string | null | undefined): PhaseInfo {
  if (!phase) return { label: 'Starting' };
  return BY_KIND[kind]?.[phase] ?? { label: humanize(phase) };
}

/** `"step 4 of 8"`, or `null` for a phase the vocabulary does not know. */
export function phaseStep(
  kind: JobKind,
  phase: string | null | undefined,
): { index: number; total: number } | null {
  if (!phase) return null;
  const order = PHASE_ORDER[kind];
  if (!order) return null;
  const index = order.indexOf(phase);
  return index < 0 ? null : { index: index + 1, total: order.length };
}
