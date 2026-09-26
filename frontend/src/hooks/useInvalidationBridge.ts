/**
 * Translates the domain keys published on the invalidation bus into the
 * cache-key prefixes declared in `services/queryKeys.ts`.
 *
 * Mounted once, by `AppLayout`. Keeping the mapping here rather than in the
 * store is what lets the store stay a plain module with no React dependency.
 */

import { useEffect } from 'react';
import { useInvalidateResources } from './useResource';
import { registerInvalidator } from '../services/invalidation';
import type { InvalidationKey } from '../services/invalidation';

/**
 * Domain key -> the cache-key first segments it makes stale.
 *
 * Every family in `services/queryKeys.ts` belongs here or in the documented
 * exclusion set this file's test pins down. The whole-backend *indexes* are the
 * easy ones to forget and the expensive ones to get wrong: `genome-index` is
 * the sole source of BOTH `/genomes` (its default, assets-only view) and
 * `/species`, so leaving it out means a finished pull or build does not change
 * either list until its five-minute stale time expires or the tab is reloaded.
 * That is invisible on a public server — server mode mounts neither the actions
 * router nor the jobs router, so nothing ever publishes on this bus — and shows
 * up only on `refgenie dash`.
 */
export const QUERY_PREFIXES: Record<InvalidationKey, readonly string[]> = {
  genomes: ['genomes', 'genome', 'genome-index', 'summary'],
  assets: [
    'assets',
    'asset',
    'asset-files',
    'asset-groups',
    'asset-group',
    'relationships',
    // A build or pull rewrites how an asset is served; `/assets/:digest` reads
    // that from the staged records.
    'staged-assets',
    'staged-asset',
  ],
  aliases: ['aliases', 'alias'],
  recipes: ['recipes', 'recipe'],
  'asset-classes': ['asset-classes', 'asset-class', 'asset-class-index'],
  servers: ['remote-servers', 'configurations', 'configuration'],
  'data-channels': ['data-channels'],
  remote: ['remote-genomes', 'remote-assets'],
  jobs: ['jobs', 'job'],
};

export function useInvalidationBridge(): void {
  const invalidate = useInvalidateResources();

  useEffect(
    () =>
      registerInvalidator((keys) => {
        const prefixes = new Set<string>();
        for (const key of keys) {
          for (const prefix of QUERY_PREFIXES[key] ?? []) prefixes.add(prefix);
        }
        for (const prefix of prefixes) {
          invalidate([prefix]);
        }
      }),
    [invalidate],
  );
}
