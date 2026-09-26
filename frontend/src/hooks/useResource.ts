/**
 * React bindings for `services/resourceCache.ts`.
 *
 * Like the cache itself this file is domain-agnostic — it is the second half of
 * the piece that could move into the shared skill. `useResource` is the single
 * read primitive every `hooks/queries/*` module is written against.
 *
 * The returned shape is intentionally the small subset of a query result this
 * app actually consumes: `data`, `error`, `isPending` and `refetch`.
 */

import { useCallback, useContext, useEffect, useMemo, useRef, useSyncExternalStore } from 'react';
import { ResourceCacheContext } from '../services/resourceCacheContext';
import { hashResourceKey } from '../services/resourceCache';
import type {
  ResourceCache,
  ResourceFetcher,
  ResourceKey,
  ResourceStatus,
  RetryPolicy,
} from '../services/resourceCache';

export interface UseResourceOptions<T, S> {
  /** False parks the read: no request, and the result stays `isPending`. */
  enabled?: boolean;
  staleTime?: number;
  retry?: RetryPolicy;
  /**
   * Transform applied outside the cache entry, so several selectors can share
   * one stored payload. MUST be a stable reference (a module-level function),
   * or it is recomputed on every render.
   */
  select?: (data: T) => S;
}

export interface UseResourceResult<S> {
  data: S | undefined;
  error: Error | undefined;
  status: ResourceStatus;
  isPending: boolean;
  isFetching: boolean;
  isError: boolean;
  isSuccess: boolean;
  dataUpdatedAt: number;
  refetch: () => Promise<void>;
}

/** The cache for the current scope. Throws outside a provider. */
export function useResourceCache(): ResourceCache {
  const cache = useContext(ResourceCacheContext);
  if (!cache) {
    throw new Error('useResource must be used inside a <ResourceCacheContext.Provider>');
  }
  return cache;
}

/** `invalidate(['genomes'])` — prefix match, see `ResourceCache.invalidate`. */
export function useInvalidateResources(): (prefix: ResourceKey) => void {
  const cache = useResourceCache();
  return useCallback((prefix: ResourceKey) => cache.invalidate(prefix), [cache]);
}

export function useResource<T, S = T>(
  key: ResourceKey,
  fetcher: ResourceFetcher<T>,
  options: UseResourceOptions<T, S> = {},
): UseResourceResult<S> {
  const cache = useResourceCache();
  const { enabled = true, select } = options;
  const hash = hashResourceKey(key);

  // The key array and the fetcher closure are rebuilt on every render, so
  // neither can be an effect dependency. The hash is the identity; the rest is
  // read through refs at the moment a request actually starts.
  const keyRef = useRef(key);
  keyRef.current = key;
  const fetcherRef = useRef(fetcher);
  fetcherRef.current = fetcher;
  const optionsRef = useRef(options);
  optionsRef.current = options;

  const stableFetcher = useCallback<ResourceFetcher<T>>((ctx) => fetcherRef.current(ctx), []);

  const subscribe = useCallback(
    (listener: () => void) => cache.subscribe(keyRef.current, hash, listener),
    [cache, hash],
  );
  const getSnapshot = useCallback(() => cache.getSnapshot<T>(hash), [cache, hash]);
  const snapshot = useSyncExternalStore(subscribe, getSnapshot, getSnapshot);

  useEffect(() => {
    if (!enabled) return;
    const { staleTime, retry } = optionsRef.current;
    void cache.ensure(keyRef.current, hash, stableFetcher, { staleTime, retry });
  }, [cache, hash, enabled, stableFetcher]);

  const refetch = useCallback(() => {
    const { staleTime, retry } = optionsRef.current;
    return cache.refetch(keyRef.current, hash, stableFetcher, { staleTime, retry });
  }, [cache, hash, stableFetcher]);

  const raw = snapshot.data;
  const data = useMemo(() => {
    if (!select) return raw as unknown as S | undefined;
    if (raw === undefined) return undefined;
    return cache.select<T, S>(hash, select, raw);
  }, [cache, hash, select, raw]);

  return useMemo(
    () => ({
      data,
      error: snapshot.error,
      status: snapshot.status,
      isPending: snapshot.isPending,
      isFetching: snapshot.isFetching,
      isError: snapshot.isError,
      isSuccess: snapshot.isSuccess,
      dataUpdatedAt: snapshot.dataUpdatedAt,
      refetch,
    }),
    [data, snapshot, refetch],
  );
}
