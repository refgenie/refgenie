/**
 * The invalidation bus.
 *
 * The job store is a plain module — it cannot hold a React hook — but a job
 * reaching `succeeded` has to refresh the same lists a mutation would. So the
 * store publishes *domain* keys here and one React-side bridge translates them
 * into TanStack query keys (`hooks/useInvalidationBridge.ts`).
 *
 * Policy, stated once: we refetch after an action, we do not edit lists
 * optimistically. refgenie mutations touch the filesystem AND the database and
 * have partial-failure and rollback paths, so an optimistically-edited list
 * would routinely disagree with the server.
 */

import type { JobKind } from './contracts';

export type InvalidationKey =
  | 'genomes'
  | 'assets'
  | 'aliases'
  | 'recipes'
  | 'asset-classes'
  | 'servers'
  | 'remote'
  | 'jobs';

export type Invalidator = (keys: readonly InvalidationKey[]) => void;

let current: Invalidator | null = null;

/** Registers the single React-side listener; returns an unregister function. */
export function registerInvalidator(fn: Invalidator): () => void {
  current = fn;
  return () => {
    if (current === fn) current = null;
  };
}

export function invalidate(keys: readonly InvalidationKey[]): void {
  current?.(keys);
}

/** What a finished job made stale, by kind. */
export const INVALIDATION_FOR_JOB_KIND: Record<JobKind, readonly InvalidationKey[]> = {
  pull: ['genomes', 'assets', 'aliases'],
  build: ['genomes', 'assets', 'aliases'],
  genome_init: ['genomes', 'aliases'],
};

/** What each curation action made stale. Fired synchronously by its handler. */
export const INVALIDATION_FOR_ACTION = {
  'asset.delete': ['assets', 'genomes'],
  'genome.delete': ['genomes', 'assets', 'aliases'],
  'alias.set': ['aliases', 'genomes'],
  'alias.remove': ['aliases', 'genomes'],
  subscribe: ['servers', 'remote'],
  unsubscribe: ['servers', 'remote'],
  'asset.set_default': ['assets'],
} as const satisfies Record<string, readonly InvalidationKey[]>;
