import { afterEach, describe, expect, it, vi } from 'vitest';
import { RETRY_DELAY_MS, createResourceCache, hashResourceKey } from './resourceCache';
import { retryUnlessClientError } from '../app/resourceCache';
import { ApiError } from './http';
import type { ResourceCache, ResourceFetcher, ResourceKey, ResourceOptions } from './resourceCache';

/** Subscribe, then ensure, the way `useResource` does. */
function read<T>(
  cache: ResourceCache,
  key: ResourceKey,
  fetcher: ResourceFetcher<T>,
  options?: ResourceOptions,
) {
  const unsubscribe = cache.subscribe(key, () => {});
  const settled = cache.ensure(key, fetcher, options);
  return { unsubscribe, settled, snapshot: () => cache.getSnapshot<T>(key) };
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

const flush = () => new Promise((resolve) => setTimeout(resolve, 0));

afterEach(() => {
  vi.useRealTimers();
});

describe('resource cache', () => {
  it('hashes keys independent of object member order and undefined members', () => {
    expect(hashResourceKey(['genomes', { pageSize: 20, page: 0, q: undefined }])).toBe(
      hashResourceKey(['genomes', { page: 0, pageSize: 20 }]),
    );
  });

  it('shares one in-flight request between readers of a key', async () => {
    const cache = createResourceCache();
    const gate = deferred<string>();
    const fetcher = vi.fn(() => gate.promise);

    const first = read(cache, ['genomes', {}], fetcher);
    read(cache, ['genomes', {}], fetcher);
    expect(first.snapshot().isPending).toBe(true);

    gate.resolve('rows');
    await first.settled;
    expect(fetcher).toHaveBeenCalledTimes(1);
    expect(first.snapshot()).toEqual({ data: 'rows', error: undefined, isPending: false });
  });

  it('serves fresh data within staleTime and refetches after it', async () => {
    vi.useFakeTimers();
    const cache = createResourceCache({ staleTime: 1000 });
    const fetcher = vi.fn(async () => 'rows');
    await read(cache, ['genomes', {}], fetcher).settled;

    await cache.ensure(['genomes', {}], fetcher);
    expect(fetcher).toHaveBeenCalledTimes(1);

    vi.advanceTimersByTime(1000);
    await cache.ensure(['genomes', {}], fetcher);
    expect(fetcher).toHaveBeenCalledTimes(2);
  });

  it('refetch ignores staleTime and drops the response it superseded', async () => {
    const cache = createResourceCache({ staleTime: 60_000 });
    const first = deferred<string>();
    const second = deferred<string>();
    const fetcher = vi.fn().mockReturnValueOnce(first.promise).mockReturnValueOnce(second.promise);

    const handle = read(cache, ['genomes', {}], fetcher);
    const refetched = cache.refetch(['genomes', {}], fetcher);
    second.resolve('new');
    first.resolve('old');
    await refetched;
    await flush();

    expect(handle.snapshot().data).toBe('new');
  });

  it('never retries a 4xx under the app policy', async () => {
    const cache = createResourceCache({ retry: retryUnlessClientError });
    const fetcher = vi.fn(async () => {
      throw new ApiError({ status: 404, detail: 'no such genome', url: '/v4/genomes/x' });
    });

    const handle = read(cache, ['genome', 'x'], fetcher);
    await handle.settled;

    expect(fetcher).toHaveBeenCalledTimes(1);
    expect(handle.snapshot().error).toBeInstanceOf(ApiError);
  });

  it('retries a 5xx twice, then fails while keeping the last good data', async () => {
    vi.useFakeTimers();
    const cache = createResourceCache({ retry: retryUnlessClientError });
    let fail = false;
    const fetcher = vi.fn(async () => {
      if (fail) throw new ApiError({ status: 503, detail: 'unavailable', url: '/v4/genomes' });
      return 'rows';
    });
    const handle = read(cache, ['genomes', {}], fetcher);
    await handle.settled;

    fail = true;
    const refetched = cache.refetch(['genomes', {}], fetcher);
    await vi.advanceTimersByTimeAsync(2 * RETRY_DELAY_MS);
    await refetched;

    expect(fetcher).toHaveBeenCalledTimes(4);
    expect(handle.snapshot().data).toBe('rows');
    expect(handle.snapshot().error).toBeInstanceOf(ApiError);
  });

  it('invalidates by key prefix, refetching only entries being read', async () => {
    const cache = createResourceCache({ staleTime: 60_000 });
    const genomes = vi.fn(async () => 'genomes');
    const other = vi.fn(async () => 'other');
    const unread = vi.fn(async () => 'unread');

    await read(cache, ['genomes', { pageSize: 20 }], genomes).settled;
    await read(cache, ['asset-files', 'genomes'], other).settled;
    const idle = read(cache, ['genomes', { pageSize: 50 }], unread);
    await idle.settled;
    idle.unsubscribe();

    cache.invalidate(['genomes']);
    await flush();
    expect(genomes).toHaveBeenCalledTimes(2);
    expect(other).toHaveBeenCalledTimes(1);
    expect(unread).toHaveBeenCalledTimes(1);

    await cache.ensure(['genomes', { pageSize: 50 }], unread);
    expect(unread).toHaveBeenCalledTimes(2);
  });

  it('aborts when the last reader leaves, without recording an error', async () => {
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
    await flush();

    expect(handle.snapshot()).toEqual({ data: undefined, error: undefined, isPending: true });
  });

  it('starts a new request when a reader remounts onto an aborted one (StrictMode)', async () => {
    const cache = createResourceCache();
    const first = deferred<string>();
    const fetcher = vi.fn().mockReturnValueOnce(first.promise).mockResolvedValueOnce('second');

    read(cache, ['genomes', {}], fetcher).unsubscribe();
    const remounted = read(cache, ['genomes', {}], fetcher);
    first.reject(new Error('aborted'));
    await remounted.settled;

    expect(fetcher).toHaveBeenCalledTimes(2);
    expect(remounted.snapshot().data).toBe('second');
  });

  it('clear drops every entry and aborts what is in flight', async () => {
    const cache = createResourceCache();
    await read(cache, ['genomes', {}], async () => 'rows').settled;
    let signal: AbortSignal | undefined;
    read(cache, ['assets', {}], ({ signal: s }) => {
      signal = s;
      return new Promise<string>(() => {});
    });

    cache.clear();

    expect(signal?.aborted).toBe(true);
    expect(cache.getSnapshot(['genomes', {}]).isPending).toBe(true);
  });
});
