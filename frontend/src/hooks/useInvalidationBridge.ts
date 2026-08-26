/**
 * Translates the domain keys published on the invalidation bus into the
 * TanStack query-key prefixes declared in `services/queryKeys.ts`.
 *
 * Mounted once, by `AppLayout`. Keeping the mapping here rather than in the
 * store is what lets the store stay a plain module with no React dependency.
 */

import { useEffect } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import { registerInvalidator } from '../services/invalidation';
import type { InvalidationKey } from '../services/invalidation';

/** Domain key -> the query-key first segments it makes stale. */
const QUERY_PREFIXES: Record<InvalidationKey, readonly string[]> = {
  genomes: ['genomes', 'genome', 'summary', 'species-summary'],
  assets: ['assets', 'asset', 'asset-files', 'asset-groups', 'asset-group', 'relationships'],
  aliases: ['aliases', 'alias'],
  recipes: ['recipes', 'recipe'],
  'asset-classes': ['asset-classes', 'asset-class'],
  servers: ['remote-servers', 'configurations', 'configuration'],
  remote: ['remote-genomes', 'remote-assets'],
  jobs: ['jobs', 'job'],
};

export function useInvalidationBridge(): void {
  const queryClient = useQueryClient();

  useEffect(
    () =>
      registerInvalidator((keys) => {
        const prefixes = new Set<string>();
        for (const key of keys) {
          for (const prefix of QUERY_PREFIXES[key] ?? []) prefixes.add(prefix);
        }
        for (const prefix of prefixes) {
          queryClient.invalidateQueries({ queryKey: [prefix] });
        }
      }),
    [queryClient],
  );
}
