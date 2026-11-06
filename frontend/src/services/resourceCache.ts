/**
 * A small keyed cache for server reads. React reads it through
 * `useSyncExternalStore` in `hooks/useResource.ts`.
 */

export type ResourceKey = readonly unknown[];

export interface ResourceSnapshot<T> {
  data: T | undefined;
  /** On a failed refresh the last good `data` is kept alongside the error. */
  error: Error | undefined;
  /** Nothing loaded or failed yet, including while a read is disabled. */
  isPending: boolean;
}

export type ResourceFetcher<T> = (ctx: { signal: AbortSignal }) => Promise<T>;

/** Return true to try again. `failures` counts the failed attempts so far. */
export type RetryPolicy = (failures: number, error: unknown) => boolean;

export interface ResourceOptions {
  /** How long a successful result is served before the next read refetches. */
  staleTime?: number;
  retry?: RetryPolicy;
}

export const NO_RETRY: RetryPolicy = () => false;
export const RETRY_DELAY_MS = 1000;

export type ResourceCache = ReturnType<typeof createResourceCache>;

/** Object members are sorted, so `{a, b}` and `{b, a}` are one entry. */
export function hashResourceKey(key: ResourceKey): string {
  return JSON.stringify(key, (_, value: unknown) =>
    value && typeof value === 'object' && !Array.isArray(value)
      ? Object.fromEntries(
          Object.keys(value)
            .sort()
            .map((name) => [name, (value as Record<string, unknown>)[name]]),
        )
      : value,
  );
}

const EMPTY: ResourceSnapshot<never> = Object.freeze({
  data: undefined,
  error: undefined,
  isPending: true,
});

interface Entry {
  key: ResourceKey;
  snapshot: ResourceSnapshot<unknown>;
  listeners: Set<() => void>;
  /** Kept so invalidation can refetch what is on screen. */
  fetcher?: ResourceFetcher<unknown>;
  options: ResourceOptions;
  controller?: AbortController;
  promise?: Promise<void>;
  /** Time of the last success; 0 means stale. */
  fetchedAt: number;
}

export function createResourceCache(defaults: ResourceOptions = {}) {
  const entries = new Map<string, Entry>();

  function entryFor(key: ResourceKey): Entry {
    const hash = hashResourceKey(key);
    let entry = entries.get(hash);
    if (!entry) {
      entry = { key, snapshot: EMPTY, listeners: new Set(), options: {}, fetchedAt: 0 };
      entries.set(hash, entry);
    }
    return entry;
  }

  function publish(entry: Entry, snapshot: ResourceSnapshot<unknown>): void {
    entry.snapshot = snapshot;
    for (const listener of entry.listeners) listener();
  }

  // Clearing the promise at once matters: a StrictMode remount must start a new
  // request, not wait on the one its own unmount just aborted.
  function cancel(entry: Entry): void {
    entry.controller?.abort();
    entry.controller = undefined;
    entry.promise = undefined;
  }

  function start(
    entry: Entry,
    fetcher: ResourceFetcher<unknown>,
    options: ResourceOptions,
  ): Promise<void> {
    cancel(entry);
    const controller = new AbortController();
    const { signal } = controller;
    const retry = options.retry ?? defaults.retry ?? NO_RETRY;
    Object.assign(entry, { fetcher, options, controller });

    const run = async () => {
      for (let failures = 0; ; failures++) {
        try {
          return await fetcher({ signal });
        } catch (error) {
          if (signal.aborted || !retry(failures, error)) throw error;
        }
        await new Promise((resolve) => setTimeout(resolve, RETRY_DELAY_MS));
        if (signal.aborted) throw signal.reason;
      }
    };

    const settle = (snapshot: ResourceSnapshot<unknown>, fetchedAt: number) => {
      // A cancelled or superseded request never touches the entry.
      if (entry.controller !== controller) return;
      entry.controller = undefined;
      entry.promise = undefined;
      entry.fetchedAt = fetchedAt;
      publish(entry, snapshot);
    };

    entry.promise = run().then(
      (data) => settle({ data, error: undefined, isPending: false }, Date.now()),
      (error) => {
        const err = error instanceof Error ? error : new Error(String(error));
        settle({ data: entry.snapshot.data, error: err, isPending: false }, 0);
      },
    );
    return entry.promise;
  }

  return {
    subscribe(key: ResourceKey, listener: () => void): () => void {
      const entry = entryFor(key);
      entry.listeners.add(listener);
      return () => {
        entry.listeners.delete(listener);
        if (entry.listeners.size === 0) cancel(entry);
      };
    },

    getSnapshot<T>(key: ResourceKey): ResourceSnapshot<T> {
      return (entries.get(hashResourceKey(key))?.snapshot ?? EMPTY) as ResourceSnapshot<T>;
    },

    /** Fetch unless a request is in flight or the data is still fresh. */
    ensure<T>(key: ResourceKey, fetcher: ResourceFetcher<T>, options: ResourceOptions = {}) {
      const entry = entryFor(key);
      entry.fetcher = fetcher;
      if (entry.promise) return entry.promise;
      const staleTime = options.staleTime ?? defaults.staleTime ?? 0;
      if (entry.fetchedAt && Date.now() - entry.fetchedAt < staleTime) return Promise.resolve();
      return start(entry, fetcher, options);
    },

    /** Fetch now, replacing any request in flight. */
    refetch<T>(key: ResourceKey, fetcher: ResourceFetcher<T>, options: ResourceOptions = {}) {
      return start(entryFor(key), fetcher, options);
    },

    /** Mark stale every entry whose key starts with `prefix`; refetch the ones on screen. */
    invalidate(prefix: ResourceKey): void {
      const wanted = hashResourceKey(prefix);
      for (const entry of entries.values()) {
        if (hashResourceKey(entry.key.slice(0, prefix.length)) !== wanted) continue;
        entry.fetchedAt = 0;
        if (entry.listeners.size > 0 && entry.fetcher) {
          void start(entry, entry.fetcher, entry.options);
        }
      }
    },

    /** Drop everything and abort every request. */
    clear(): void {
      for (const entry of entries.values()) {
        cancel(entry);
        publish(entry, EMPTY);
      }
      entries.clear();
    },
  };
}
