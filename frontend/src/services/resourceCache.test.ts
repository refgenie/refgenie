/**
 * The generic cache, exercised without React.
 *
 * Every behaviour here is one the app depended on TanStack Query for, so this
 * file is the contract that replaced it: shared entries, in-flight dedup, stale
 * time, retry policy, prefix invalidation, cancellation and memoized selects.
 */

import { describe, expect, it, vi } from 'vitest';
import { NO_RETRY, createResourceCache, hashResourceKey } from './resourceCache';
import { retryUnlessClientError } from '../app/resourceCache';
import { ApiError } from './http';
import type { ResourceCache, ResourceFetcher, ResourceKey } from './resourceCache';

/** Subscribe + ensure in one call, the way `useResource` does it. */
function read<T>(
  cache: ResourceCache,
  key: ResourceKey,
  fetcher: ResourceFetcher<T>,
  options?: Parameters<ResourceCache['ensure']>[3],
) {
  const hash = hashResourceKey(key);
  const unsubscribe = cache.subscribe(key, hash, () => {});
  return {
    hash,
    unsubscribe,
    settled: cache.ensure(key, hash, fetcher, options),
    snapshot: () => cache.getSnapshot<T>(hash),
  };
}

function deferred<T>() {
  let resolve!: (value: T) => void;
  let reject!: (error: unknown) => void;
  const promise = new Promise<T>((res, rej) => {
    resolve = res;
    reject = rej;
  });
  return { promise, resolve, reject };
}

const noDelay = { retryDelay: () => 0 };

describe('hashResourceKey', () => {
  it('is insensitive to object key order', () => {
    expect(hashResourceKey(['genomes', { limit: 20, offset: 0 }])).toBe(
      hashResourceKey(['genomes', { offset: 0, limit: 20 }]),
    );
  });

  it('treats an undefined member as an absent one', () => {
    expect(hashResourceKey(['genomes', { limit: 20, q: undefined }])).toBe(
      hashResourceKey(['genomes', { limit: 20 }]),
    );
  });

  it('separates different keys', () => {
    expect(hashResourceKey(['genomes', { limit: 20 }])).not.toBe(
      hashResourceKey(['genomes', { limit: 50 }]),
    );
  });
});

describe('resource cache', () => {
  it('serves two readers of one key from a single request', async () => {
    const cache = createResourceCache();
    const fetcher = vi.fn(async () => 'rows');

    const first = read(cache, ['genomes', {}], fetcher, { staleTime: 60_000 });
    await first.settled;
    const second = read(cache, ['genomes', {}], fetcher, { staleTime: 60_000 });
    await second.settled;

    expect(fetcher).toHaveBeenCalledTimes(1);
    expect(first.snapshot().data).toBe('rows');
    expect(second.snapshot()).toBe(first.snapshot());
  });

  it('deduplicates concurrent requests for one key', async () => {
    const cache = createResourceCache();
    const gate = deferred<string>();
    const fetcher = vi.fn(() => gate.promise);

    const first = read(cache, ['genomes', {}], fetcher);
    const second = read(cache, ['genomes', {}], fetcher);
    expect(fetcher).toHaveBeenCalledTimes(1);
    expect(first.snapshot().isFetching).toBe(true);

    gate.resolve('rows');
    await Promise.all([first.settled, second.settled]);

    expect(fetcher).toHaveBeenCalledTimes(1);
    expect(second.snapshot().data).toBe('rows');
  });

  it('reports pending, then success', async () => {
    const cache = createResourceCache();
    const gate = deferred<string>();
    const handle = read(cache, ['genome', 'abc'], () => gate.promise);

    expect(handle.snapshot().status).toBe('pending');
    expect(handle.snapshot().isPending).toBe(true);

    gate.resolve('one');
    await handle.settled;

    expect(handle.snapshot().status).toBe('success');
    expect(handle.snapshot().isPending).toBe(false);
    expect(handle.snapshot().isFetching).toBe(false);
  });

  it('serves a fresh entry without refetching, and refetches once stale', async () => {
    const cache = createResourceCache();
    const fetcher = vi.fn(async () => 'rows');

    const first = read(cache, ['genomes', {}], fetcher, { staleTime: 10 });
    await first.settled;
    expect(fetcher).toHaveBeenCalledTimes(1);

    // Inside the stale window: cache hit, no request.
    await cache.ensure(['genomes', {}], first.hash, fetcher, { staleTime: 10 });
    expect(fetcher).toHaveBeenCalledTimes(1);

    await new Promise((resolve) => setTimeout(resolve, 15));
    await cache.ensure(['genomes', {}], first.hash, fetcher, { staleTime: 10 });
    expect(fetcher).toHaveBeenCalledTimes(2);
  });

  it('honours a long stale time, which is what makes the index one request', async () => {
    const cache = createResourceCache({ staleTime: 5 * 60 * 1000 });
    const fetcher = vi.fn(async () => [1, 2, 3]);

    const first = read(cache, ['genome-index'], fetcher);
    await first.settled;
    first.unsubscribe();

    const second = read(cache, ['genome-index'], fetcher);
    await second.settled;

    expect(fetcher).toHaveBeenCalledTimes(1);
    expect(second.snapshot().data).toEqual([1, 2, 3]);
  });

  it('refetch ignores the stale window', async () => {
    const cache = createResourceCache({ staleTime: 60_000 });
    const fetcher = vi.fn(async () => 'rows');
    const handle = read(cache, ['genomes', {}], fetcher);
    await handle.settled;

    await cache.refetch(['genomes', {}], handle.hash, fetcher);
    expect(fetcher).toHaveBeenCalledTimes(2);
  });
});

describe('retry policy', () => {
  it('never retries a 4xx', async () => {
    const cache = createResourceCache({ retry: retryUnlessClientError, ...noDelay });
    const fetcher = vi.fn(async () => {
      throw new ApiError({ status: 404, detail: 'no such genome', url: '/v4/genomes/x' });
    });

    const handle = read(cache, ['genome', 'x'], fetcher);
    await handle.settled;

    expect(fetcher).toHaveBeenCalledTimes(1);
    expect(handle.snapshot().status).toBe('error');
    expect(handle.snapshot().error).toBeInstanceOf(ApiError);
  });

  it('retries a 5xx twice, then gives up', async () => {
    const cache = createResourceCache({ retry: retryUnlessClientError, ...noDelay });
    const fetcher = vi.fn(async () => {
      throw new ApiError({ status: 503, detail: 'unavailable', url: '/v4/genomes' });
    });

    const handle = read(cache, ['genomes', {}], fetcher);
    await handle.settled;

    expect(fetcher).toHaveBeenCalledTimes(3);
    expect(handle.snapshot().isError).toBe(true);
  });

  it('retries a transport failure, which carries status 0', async () => {
    const cache = createResourceCache({ retry: retryUnlessClientError, ...noDelay });
    let attempts = 0;
    const fetcher = vi.fn(async () => {
      attempts += 1;
      if (attempts < 3) {
        throw new ApiError({ status: 0, detail: 'offline', url: '/v4/genomes' });
      }
      return 'rows';
    });

    const handle = read(cache, ['genomes', {}], fetcher);
    await handle.settled;

    expect(fetcher).toHaveBeenCalledTimes(3);
    expect(handle.snapshot().data).toBe('rows');
  });

  it('NO_RETRY gives up immediately', async () => {
    const cache = createResourceCache({ ...noDelay });
    const fetcher = vi.fn(async () => {
      throw new ApiError({ status: 500, detail: 'boom', url: '/v4/summary' });
    });

    const handle = read(cache, ['summary'], fetcher, { retry: NO_RETRY });
    await handle.settled;

    expect(fetcher).toHaveBeenCalledTimes(1);
  });

  it('keeps the last good rows when a refresh fails', async () => {
    const cache = createResourceCache({ ...noDelay });
    let fail = false;
    const fetcher = async () => {
      if (fail) throw new ApiError({ status: 500, detail: 'boom', url: '/v4/genomes' });
      return 'rows';
    };

    const handle = read(cache, ['genomes', {}], fetcher, { retry: NO_RETRY });
    await handle.settled;
    fail = true;
    await cache.refetch(['genomes', {}], handle.hash, fetcher, { retry: NO_RETRY });

    expect(handle.snapshot().data).toBe('rows');
    expect(handle.snapshot().status).toBe('error');
  });
});

describe('invalidation', () => {
  it('matches by key prefix and refreshes what is being read', async () => {
    const cache = createResourceCache({ staleTime: 60_000 });
    const genomes = vi.fn(async () => 'genomes');
    const assets = vi.fn(async () => 'assets');

    const page1 = read(cache, ['genomes', { limit: 20 }], genomes);
    const page2 = read(cache, ['genomes', { limit: 50 }], genomes);
    const assetList = read(cache, ['assets', {}], assets);
    await Promise.all([page1.settled, page2.settled, assetList.settled]);
    expect(genomes).toHaveBeenCalledTimes(2);

    cache.invalidate(['genomes']);
    await Promise.resolve();
    await new Promise((resolve) => setTimeout(resolve, 0));

    expect(genomes).toHaveBeenCalledTimes(4);
    expect(assets).toHaveBeenCalledTimes(1);
  });

  it('does not match a key that merely shares a segment further along', async () => {
    const cache = createResourceCache({ staleTime: 60_000 });
    const fetcher = vi.fn(async () => 'rows');
    const handle = read(cache, ['asset-files', 'genomes'], fetcher);
    await handle.settled;

    cache.invalidate(['genomes']);
    await new Promise((resolve) => setTimeout(resolve, 0));

    expect(fetcher).toHaveBeenCalledTimes(1);
  });

  it('defers an unread entry to its next read rather than refetching it now', async () => {
    const cache = createResourceCache({ staleTime: 60_000 });
    const fetcher = vi.fn(async () => 'rows');
    const handle = read(cache, ['genomes', {}], fetcher);
    await handle.settled;
    handle.unsubscribe();

    cache.invalidate(['genomes']);
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(fetcher).toHaveBeenCalledTimes(1);

    const again = read(cache, ['genomes', {}], fetcher);
    await again.settled;
    expect(fetcher).toHaveBeenCalledTimes(2);
  });
});

describe('cancellation', () => {
  it('passes an AbortSignal to the fetcher', async () => {
    const cache = createResourceCache();
    let seen: AbortSignal | undefined;
    const handle = read(cache, ['genomes', {}], async ({ signal }) => {
      seen = signal;
      return 'rows';
    });
    await handle.settled;

    expect(seen).toBeInstanceOf(AbortSignal);
    expect(seen?.aborted).toBe(false);
  });

  it('aborts when the last reader goes away, and leaves no error behind', async () => {
    const cache = createResourceCache();
    const gate = deferred<string>();
    let signal: AbortSignal | undefined;
    const handle = read(cache, ['genomes', {}], ({ signal: s }) => {
      signal = s;
      return gate.promise;
    });

    handle.unsubscribe();
    expect(signal?.aborted).toBe(true);

    gate.reject(new Error('aborted'));
    await handle.settled;

    expect(handle.snapshot().status).toBe('pending');
    expect(handle.snapshot().error).toBeUndefined();
    expect(handle.snapshot().isFetching).toBe(false);
  });

  it('starts a fresh request when a reader remounts onto an aborted one', async () => {
    // StrictMode's mount/unmount/remount, which lands the second `ensure`
    // inside the window where the first request is aborted but its rejection
    // has not been delivered yet. Deduping onto it strands the entry at
    // `pending` with nothing in flight — a skeleton that never resolves.
    const cache = createResourceCache();
    const first = deferred<string>();
    let call = 0;
    const fetcher = () => (call++ === 0 ? first.promise : Promise.resolve('second'));

    const mounted = read(cache, ['genomes', {}], fetcher);
    mounted.unsubscribe();

    const remounted = read(cache, ['genomes', {}], fetcher);
    first.reject(new Error('aborted'));
    await remounted.settled;

    expect(call).toBe(2);
    expect(remounted.snapshot().status).toBe('success');
    expect(remounted.snapshot().data).toBe('second');
  });

  it('drops a superseded response in favour of the refetch that replaced it', async () => {
    const cache = createResourceCache();
    const first = deferred<string>();
    const second = deferred<string>();
    let call = 0;
    const fetcher = () => (call++ === 0 ? first.promise : second.promise);

    const handle = read(cache, ['genomes', {}], fetcher);
    const refetched = cache.refetch(['genomes', {}], handle.hash, fetcher);

    second.resolve('new');
    first.resolve('old');
    await refetched;

    expect(handle.snapshot().data).toBe('new');
  });
});

describe('select', () => {
  it('computes each selector once and shares the result', async () => {
    const cache = createResourceCache();
    const toUpper = vi.fn((rows: string[]) => rows.map((r) => r.toUpperCase()));
    const toCount = vi.fn((rows: string[]) => rows.length);

    const handle = read(cache, ['genome-index'], async () => ['a', 'b']);
    await handle.settled;
    const rows = handle.snapshot().data as string[];

    const upperOnce = cache.select(handle.hash, toUpper, rows);
    const upperTwice = cache.select(handle.hash, toUpper, rows);
    const count = cache.select(handle.hash, toCount, rows);

    expect(toUpper).toHaveBeenCalledTimes(1);
    expect(toCount).toHaveBeenCalledTimes(1);
    expect(upperTwice).toBe(upperOnce);
    expect(count).toBe(2);
  });

  it('recomputes when the underlying rows are replaced', async () => {
    const cache = createResourceCache();
    const select = vi.fn((rows: string[]) => rows.length);
    const handle = read(cache, ['genome-index'], async () => ['a']);
    await handle.settled;

    cache.select(handle.hash, select, handle.snapshot().data as string[]);
    cache.select(handle.hash, select, ['a', 'b']);

    expect(select).toHaveBeenCalledTimes(2);
  });
});

describe('window focus', () => {
  it('does not refetch when the window regains focus', async () => {
    const cache = createResourceCache({ staleTime: 60_000 });
    const fetcher = vi.fn(async () => 'rows');
    const handle = read(cache, ['genomes', {}], fetcher);
    await handle.settled;

    window.dispatchEvent(new Event('focus'));
    window.dispatchEvent(new Event('visibilitychange'));
    document.dispatchEvent(new Event('visibilitychange'));
    await new Promise((resolve) => setTimeout(resolve, 0));

    expect(fetcher).toHaveBeenCalledTimes(1);
  });
});

describe('clear', () => {
  it('drops every entry and cancels what is in flight', async () => {
    const cache = createResourceCache({ staleTime: 60_000 });
    const gate = deferred<string>();
    let signal: AbortSignal | undefined;
    const settled = read(cache, ['genomes', {}], async () => 'rows');
    await settled.settled;
    const inFlight = read(cache, ['assets', {}], ({ signal: s }) => {
      signal = s;
      return gate.promise;
    });

    cache.clear();

    expect(signal?.aborted).toBe(true);
    expect(cache.peek(settled.hash)).toBeUndefined();
    expect(cache.getSnapshot(inFlight.hash).data).toBeUndefined();
  });
});
