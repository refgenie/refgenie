import { QueryClient } from '@tanstack/react-query';
import { ApiError } from '../services/http';

/**
 * Explicitly configured, never a bare `new QueryClient()`.
 *
 * A 4xx is a client mistake and is never worth retrying; a 5xx gets two more
 * attempts.
 */
export function createQueryClient(): QueryClient {
  return new QueryClient({
    defaultOptions: {
      queries: {
        staleTime: 30_000,
        refetchOnWindowFocus: false,
        retry: (failureCount, error) => {
          if (error instanceof ApiError && error.status >= 400 && error.status < 500) {
            return false;
          }
          return failureCount < 2;
        },
      },
    },
  });
}
