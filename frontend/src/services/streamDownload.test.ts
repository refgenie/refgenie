// @vitest-environment node
import { readFileSync } from 'node:fs';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { downloadUrl, servePort } from './streamDownload';
import type { PortLike } from './streamDownload';
import { BLOB_FALLBACK_MAX, chooseSaveMethod } from './saveStream';

type Handler = (event: unknown) => void;

/** Run public/download-sw.js against a fake worker global. */
function loadWorker(scope = 'http://localhost/app/') {
  const handlers: Record<string, Handler> = {};
  const self = {
    addEventListener: (type: string, fn: Handler) => (handlers[type] = fn),
    registration: { scope },
    skipWaiting: () => {},
    clients: { claim: async () => {} },
  };
  const code = readFileSync(
    new URL('../../public/download-sw.js', import.meta.url),
    'utf8',
  );
  new Function('self', code)(self);
  return {
    register(id: string, filename: string, port: MessagePort) {
      handlers.message({ data: { type: 'register', id, filename }, ports: [port] });
    },
    fetch(url: string): Promise<Response> | undefined {
      let response: Promise<Response> | undefined;
      handlers.fetch({
        request: { url },
        respondWith: (r: Response | Promise<Response>) =>
          (response = Promise.resolve(r)),
      });
      return response;
    },
  };
}

async function* bytes(parts: number[][], failAfter?: number) {
  for (let i = 0; i < parts.length; i += 1) {
    if (failAfter === i) throw new Error('store went away');
    yield Uint8Array.from(parts[i]);
  }
}

function start(
  source: AsyncIterable<Uint8Array>,
  opts: { signal?: AbortSignal; scope?: string } = {},
) {
  const worker = loadWorker(opts.scope);
  const channel = new MessageChannel();
  const server = servePort(channel.port1 as unknown as PortLike, source, {
    signal: opts.signal,
    batchBytes: 4,
  });
  worker.register('abc', 'hg38.fa.gz', channel.port2);
  return { worker, server };
}

describe('download worker and page', () => {
  it('streams every byte, in order, as an attachment', async () => {
    const { worker, server } = start(bytes([[1, 2], [3], [4, 5, 6], [7]]));
    await server.registered;
    const response = await worker.fetch(
      'http://localhost/app/rg-download/abc/hg38.fa.gz',
    );
    expect(response).toBeDefined();
    await server.started;
    expect(response!.headers.get('Content-Disposition')).toMatch(
      /^attachment; filename="hg38.fa.gz"/,
    );
    const body = new Uint8Array(await response!.arrayBuffer());
    expect([...body]).toEqual([1, 2, 3, 4, 5, 6, 7]);
    await expect(server.done).resolves.toBeUndefined();
  });

  it('leaves other requests alone', () => {
    const worker = loadWorker();
    expect(worker.fetch('http://localhost/app/genomes/abc')).toBeUndefined();
  });

  it('answers an unknown or used download id with 404', async () => {
    const { worker, server } = start(bytes([[1]]));
    await server.registered;
    await worker
      .fetch('http://localhost/app/rg-download/abc/x')!
      .then((r) => r.arrayBuffer());
    const again = await worker.fetch('http://localhost/app/rg-download/abc/x');
    expect(again!.status).toBe(404);
  });

  it('errors the response stream when the source fails, so no complete-looking file', async () => {
    const { worker, server } = start(bytes([[1, 2, 3, 4], [5, 6, 7, 8], [9]], 2));
    await server.registered;
    const response = await worker.fetch('http://localhost/app/rg-download/abc/x');
    await expect(response!.arrayBuffer()).rejects.toThrow();
    await expect(server.done).rejects.toThrow('store went away');
  });

  it('errors the response stream when the page cancels', async () => {
    const abort = new AbortController();
    async function* slow() {
      yield Uint8Array.from([1, 2, 3, 4]);
      abort.abort();
      await new Promise((resolve) => setTimeout(resolve, 10));
      yield Uint8Array.from([5]);
    }
    const { worker, server } = start(slow(), { signal: abort.signal });
    await server.registered;
    const response = await worker.fetch('http://localhost/app/rg-download/abc/x');
    await expect(response!.arrayBuffer()).rejects.toThrow();
    await expect(server.done).rejects.toMatchObject({ name: 'AbortError' });
  });

  it('stops the source when the browser cancels the download', async () => {
    let produced = 0;
    async function* endless() {
      for (;;) {
        produced += 1;
        yield new Uint8Array(4);
      }
    }
    const { worker, server } = start(endless());
    await server.registered;
    const response = await worker.fetch('http://localhost/app/rg-download/abc/x');
    const reader = response!.body!.getReader();
    await reader.read();
    await reader.cancel();
    await expect(server.done).rejects.toMatchObject({ name: 'AbortError' });
    const after = produced;
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(produced).toBe(after);
  });

  it('only pulls from the source when the browser reads (backpressure)', async () => {
    let produced = 0;
    async function* counted() {
      for (let i = 0; i < 1000; i += 1) {
        produced += 1;
        yield new Uint8Array(4);
      }
    }
    const { worker, server } = start(counted());
    await server.registered;
    const response = await worker.fetch('http://localhost/app/rg-download/abc/x');
    await new Promise((resolve) => setTimeout(resolve, 50));
    // A batch is 4 bytes here, one chunk; the worker buffers only a couple.
    expect(produced).toBeLessThan(10);
    await response!.body!.cancel();
  });
});

const SCOPE = 'http://localhost/app/rg-download/';

describe('download worker at its own narrow scope', () => {
  afterEach(() => {
    vi.useRealTimers();
  });

  it('builds download URLs under the scope', () => {
    expect(downloadUrl(SCOPE, 'abc', 'hg 38.fa.gz')).toBe(
      'http://localhost/app/rg-download/abc/hg%2038.fa.gz',
    );
  });

  it('streams a download under the scope and leaves the rest alone', async () => {
    const { worker, server } = start(bytes([[1, 2], [3]]), { scope: SCOPE });
    await server.registered;
    expect(worker.fetch('http://localhost/app/genomes/abc')).toBeUndefined();
    const response = await worker.fetch(downloadUrl(SCOPE, 'abc', 'x.fa.gz'));
    expect([...new Uint8Array(await response!.arrayBuffer())]).toEqual([1, 2, 3]);
    await expect(server.done).resolves.toBeUndefined();
  });

  it('serves the newest registration after the page registers again', async () => {
    const { worker, server } = start(bytes([[1, 2], [3]]), { scope: SCOPE });
    await server.registered;
    const again = new MessageChannel();
    const reRegistered = server.adopt(again.port1 as unknown as PortLike);
    worker.register('abc', 'x.fa.gz', again.port2);
    await reRegistered;
    const response = await worker.fetch(downloadUrl(SCOPE, 'abc', 'x.fa.gz'));
    await server.started;
    expect([...new Uint8Array(await response!.arrayBuffer())]).toEqual([1, 2, 3]);
    await expect(server.done).resolves.toBeUndefined();
  });

  it('drops a registration the page gives up on before it is fetched', async () => {
    const { worker, server } = start(bytes([[1]]), { scope: SCOPE });
    await server.registered;
    server.stop(new Error('Not started'));
    await expect(server.done).rejects.toThrow('Not started');
    await new Promise((resolve) => setTimeout(resolve, 10));
    const response = await worker.fetch(downloadUrl(SCOPE, 'abc', 'x'));
    expect(response!.status).toBe(404);
  });

  it('drops a registration that is never fetched', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout'] });
    const { worker, server } = start(bytes([[1]]), { scope: SCOPE });
    await server.registered;
    vi.advanceTimersByTime(60_000);
    const response = await worker.fetch(downloadUrl(SCOPE, 'abc', 'x'));
    expect(response!.status).toBe(404);
    server.stop(new Error('done'));
  });
});

describe('servePort.stop', () => {
  it('removes its abort listener and stops the source', async () => {
    const abort = new AbortController();
    const remove = vi.spyOn(abort.signal, 'removeEventListener');
    let returned = false;
    const source: AsyncIterable<Uint8Array> = {
      [Symbol.asyncIterator]: () => ({
        next: async () => ({ done: false, value: new Uint8Array(1) }),
        return: async () => {
          returned = true;
          return { done: true, value: undefined };
        },
      }),
    };
    const channel = new MessageChannel();
    const server = servePort(channel.port1 as unknown as PortLike, source, {
      signal: abort.signal,
    });
    server.stop(new Error('Not started'));
    expect(remove).toHaveBeenCalledWith('abort', expect.any(Function));
    abort.abort();
    await expect(server.done).rejects.toThrow('Not started');
    expect(returned).toBe(true);
    channel.port2.close();
  });
});

describe('chooseSaveMethod', () => {
  it('prefers the file picker, then the stream, then a small Blob', () => {
    expect(
      chooseSaveMethod({ hasFilePicker: true, canStream: true, bytes: 1e10 }),
    ).toBe('file-picker');
    expect(
      chooseSaveMethod({ hasFilePicker: false, canStream: true, bytes: 1e10 }),
    ).toBe('stream');
    expect(
      chooseSaveMethod({ hasFilePicker: false, canStream: false, bytes: 1000 }),
    ).toBe('blob');
    expect(
      chooseSaveMethod({
        hasFilePicker: false,
        canStream: false,
        bytes: BLOB_FALLBACK_MAX + 1,
      }),
    ).toBeNull();
  });
});
