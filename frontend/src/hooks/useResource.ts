import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useSyncExternalStore } from 'react';
import { hashResourceKey } from '../services/resourceCache';
import type {
  ResourceCache,
  ResourceFetcher,
  ResourceKey,
  ResourceOptions,
} from '../services/resourceCache';

export const ResourceCacheContext = createContext<ResourceCache | null>(null);

export interface UseResourceOptions extends ResourceOptions {
  /** False skips the request; the result stays pending. */
  enabled?: boolean;
}

export function useResourceCache(): ResourceCache {
  const cache = useContext(ResourceCacheContext);
  if (!cache) throw new Error('useResource needs a ResourceCacheContext provider');
  return cache;
}

/** Returns a function that marks stale every key starting with one of `families`. */
export function useInvalidate(): (families: readonly string[]) => void {
  const cache = useResourceCache();
  return useCallback(
    (families: readonly string[]) => {
      for (const family of families) cache.invalidate([family]);
    },
    [cache],
  );
}

export function useResource<T>(
  key: ResourceKey,
  fetcher: ResourceFetcher<T>,
  { enabled = true, staleTime, retry }: UseResourceOptions = {},
) {
  const cache = useResourceCache();
  // Callers rebuild the key array and fetcher every render; the hash is the identity.
  const hash = hashResourceKey(key);
  const stableKey = useMemo(() => JSON.parse(hash) as ResourceKey, [hash]);
  const latest = useRef({ fetcher, options: { staleTime, retry } });
  latest.current = { fetcher, options: { staleTime, retry } };
  const read = useCallback<ResourceFetcher<T>>((ctx) => latest.current.fetcher(ctx), []);

  const snapshot = useSyncExternalStore(
    useCallback((listener: () => void) => cache.subscribe(stableKey, listener), [cache, stableKey]),
    () => cache.getSnapshot<T>(stableKey),
  );

  useEffect(() => {
    if (enabled) void cache.ensure(stableKey, read, latest.current.options);
  }, [cache, stableKey, enabled, read]);

  const refetch = useCallback(
    () => cache.refetch(stableKey, read, latest.current.options),
    [cache, stableKey, read],
  );

  return useMemo(() => ({ ...snapshot, refetch }), [snapshot, refetch]);
}
