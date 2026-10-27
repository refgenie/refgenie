/**
 * The app's cache policy, in one place — the successor to `createQueryClient()`.
 *
 * Explicitly configured, never a bare `createResourceCache()`. A 4xx is a
 * client mistake and is never worth retrying; a 5xx or a transport failure gets
 * two more attempts.
 *
 * There is no window-focus refetching anywhere in this app: the cache has no
 * focus listener at all, so a stale entry refreshes when something reads it,
 * when a `staleTime` lapses, or when the invalidation bus says so — never
 * because the user alt-tabbed.
 */

import { createResourceCache } from '../services/resourceCache';
import { ApiError } from '../services/http';
import type { ResourceCache, RetryPolicy } from '../services/resourceCache';

/** Every read is served from cache for this long before a refresh. */
export const DEFAULT_STALE_MS = 30_000;

export const retryUnlessClientError: RetryPolicy = (failureCount, error) => {
  if (error instanceof ApiError && error.status >= 400 && error.status < 500) {
    return false;
  }
  return failureCount < 2;
};

export function createAppResourceCache(): ResourceCache {
  return createResourceCache({
    staleTime: DEFAULT_STALE_MS,
    retry: retryUnlessClientError,
  });
}
