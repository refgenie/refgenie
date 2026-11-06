import { createResourceCache } from '../services/resourceCache';
import { ApiError } from '../services/http';
import type { RetryPolicy } from '../services/resourceCache';

/** A 4xx is the client's mistake and will fail again; anything else gets two retries. */
export const retryUnlessClientError: RetryPolicy = (failures, error) =>
  !(error instanceof ApiError && error.status >= 400 && error.status < 500) && failures < 2;

export const createAppResourceCache = () =>
  createResourceCache({ staleTime: 30_000, retry: retryUnlessClientError });

/** The root cache, exported so the job store can invalidate it. */
export const appCache = createAppResourceCache();
