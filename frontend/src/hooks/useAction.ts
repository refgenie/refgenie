import { useCallback, useState } from 'react';
import { ApiError } from '../services/http';

/**
 * Runs a mutation with a busy flag. A failure becomes `error` (the server's
 * `detail` when it sent one), or goes to `onError`, for a toast, when given.
 * `run` resolves true when `fn` succeeded.
 */
export function useAction(onError?: (message: string) => void, fallback?: string) {
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const run = useCallback(
    async (fn: () => Promise<unknown>): Promise<boolean> => {
      setBusy(true);
      setError(null);
      try {
        await fn();
        return true;
      } catch (caught) {
        const message =
          caught instanceof ApiError ? caught.detail : (fallback ?? String(caught));
        if (onError) onError(message);
        else setError(message);
        return false;
      } finally {
        setBusy(false);
      }
    },
    [onError, fallback],
  );

  return { busy, error, setError, run };
}
