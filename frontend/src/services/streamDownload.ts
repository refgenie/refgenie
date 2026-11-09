/**
 * Stream a download to disk through a service worker, for browsers without
 * the File System Access API (Firefox, Safari).
 *
 * The page registers the download with `public/download-sw.js` over a
 * MessageChannel, then opens a unique URL in a hidden iframe. The worker
 * answers that URL with a streamed attachment and asks the page for one batch
 * of bytes each time the browser wants more, so only a few batches are ever
 * in memory. Keep the message protocol in step with the worker.
 */

/** Bytes sent per message to the worker. */
export const BATCH_BYTES = 1 << 20;
/**
 * The worker's scope, relative to the page base. Only download URLs live
 * here, so the worker controls nothing else on the site.
 */
export const DOWNLOAD_PATH = 'rg-download/';
const WORKER_FILE = 'download-sw.js';
const START_TIMEOUT_MS = 15_000;
const PING_MS = 10_000;

/** What the worker sends on the channel. */
export type WorkerMessage =
  { type: 'registered' } | { type: 'started' } | { type: 'pull' } | { type: 'cancel' };

/** What the page sends on the channel. */
export type PageMessage =
  | { type: 'chunk'; chunk: Uint8Array }
  | { type: 'end' }
  | { type: 'error'; message: string };

/** The page's end of the channel: a MessagePort, or a fake in tests. */
export interface PortLike {
  postMessage(message: PageMessage, transfer?: Transferable[]): void;
  onmessage: ((event: { data: WorkerMessage }) => void) | null;
  close(): void;
}

export interface PortServer {
  /** The worker has the download and its URL can be opened. */
  registered: Promise<void>;
  /**
   * Move to a fresh port, for registering the download again; resolves when
   * the worker acknowledges it there. The old port is closed.
   */
  adopt(port: PortLike): Promise<void>;
  /** Give up: tell the worker, stop the source, and reject `done`. */
  stop(error: unknown): void;
  /** The browser has asked the worker for the URL: the download has begun. */
  started: Promise<void>;
  /** All bytes sent and the stream ended; rejects on failure or cancel. */
  done: Promise<void>;
}

function abortError(message: string): DOMException {
  return new DOMException(message, 'AbortError');
}

/**
 * Answer the worker's pulls from `source`, batching chunks into messages of
 * about `batchBytes`. On a source error or an abort the worker is told to
 * error the stream, so the browser marks the download as failed.
 */
export function servePort(
  port: PortLike,
  source: AsyncIterable<Uint8Array>,
  {
    signal,
    batchBytes = BATCH_BYTES,
  }: { signal?: AbortSignal; batchBytes?: number } = {},
): PortServer {
  const iterator = source[Symbol.asyncIterator]();
  let current = port;
  let onRegistered!: () => void;
  let onStarted!: () => void;
  let resolveDone!: () => void;
  let rejectDone!: (error: unknown) => void;
  const registered = new Promise<void>((resolve) => (onRegistered = resolve));
  const started = new Promise<void>((resolve) => (onStarted = resolve));
  const done = new Promise<void>((resolve, reject) => {
    resolveDone = resolve;
    rejectDone = reject;
  });
  done.catch(() => {});

  let finished = false;
  let sourceDone = false;

  const finish = (error?: unknown, tellWorker = true) => {
    if (finished) return;
    finished = true;
    signal?.removeEventListener('abort', onAbort);
    if (error === undefined) {
      current.postMessage({ type: 'end' });
      resolveDone();
    } else {
      if (tellWorker) {
        const message = error instanceof Error ? error.message : String(error);
        current.postMessage({ type: 'error', message });
      }
      iterator.return?.().catch(() => {});
      rejectDone(error);
    }
    current.onmessage = null;
    current.close();
  };

  const onAbort = () => finish(signal?.reason ?? abortError('Cancelled'));
  if (signal?.aborted) {
    onAbort();
  } else {
    signal?.addEventListener('abort', onAbort, { once: true });
  }

  const sendBatch = async () => {
    const parts: Uint8Array[] = [];
    let size = 0;
    while (!sourceDone && size < batchBytes) {
      const next = await iterator.next();
      if (finished) return;
      if (next.done) {
        sourceDone = true;
        break;
      }
      parts.push(next.value);
      size += next.value.length;
    }
    if (size === 0) {
      finish();
      return;
    }
    // Always a fresh copy: the buffer is transferred, and a source may reuse its own.
    const batch = concat(parts, size);
    current.postMessage({ type: 'chunk', chunk: batch }, [batch.buffer]);
  };

  const listen = (target: PortLike) => {
    target.onmessage = (event) => {
      if (finished || target !== current) return;
      const msg = event.data;
      if (msg.type === 'registered') onRegistered();
      else if (msg.type === 'started') onStarted();
      else if (msg.type === 'pull') sendBatch().catch((error: unknown) => finish(error));
      else if (msg.type === 'cancel')
        finish(abortError('The download was cancelled'), false);
    };
  };
  listen(port);

  const adopt = (next: PortLike): Promise<void> => {
    if (finished) {
      next.close();
      return done;
    }
    current.onmessage = null;
    current.close();
    current = next;
    const again = new Promise<void>((resolve) => (onRegistered = resolve));
    listen(next);
    return again;
  };

  return { registered, started, done, adopt, stop: (error) => finish(error) };
}

function concat(parts: Uint8Array[], size: number): Uint8Array {
  const out = new Uint8Array(size);
  let offset = 0;
  for (const part of parts) {
    out.set(part, offset);
    offset += part.length;
  }
  return out;
}

/** True if this page can stream a download through a service worker. */
export function canStreamDownload(): boolean {
  return (
    typeof window !== 'undefined' &&
    window.isSecureContext &&
    typeof navigator !== 'undefined' &&
    'serviceWorker' in navigator &&
    typeof MessageChannel !== 'undefined'
  );
}

let workerPromise: Promise<ServiceWorkerRegistration> | null = null;

/** Register the download worker once and wait until it is active. */
function downloadWorker(): Promise<ServiceWorkerRegistration> {
  workerPromise ??= (async () => {
    // Relative to <base href>, so it works behind a path prefix too.
    const url = new URL(WORKER_FILE, document.baseURI);
    const registration = await navigator.serviceWorker.register(url.href, {
      scope: new URL(DOWNLOAD_PATH, document.baseURI).href,
    });
    const worker = registration.installing ?? registration.waiting;
    if (!registration.active && worker) {
      await new Promise<void>((resolve, reject) => {
        const check = () => {
          if (worker.state === 'activated') resolve();
          else if (worker.state === 'redundant')
            reject(new Error('The download helper did not start.'));
        };
        worker.addEventListener('statechange', check);
        check();
      });
    }
    return registration;
  })();
  workerPromise.catch(() => {
    workerPromise = null;
  });
  return workerPromise;
}

/** The URL, under the worker's scope, that the worker answers for download `id`. */
export function downloadUrl(scope: string, id: string, filename: string): string {
  return new URL(`${id}/${encodeURIComponent(filename)}`, scope).href;
}

function timeout(ms: number, message: string): Promise<never> {
  return new Promise((_, reject) => setTimeout(() => reject(new Error(message)), ms));
}

/**
 * Unregister a download worker left at the old site-wide scope, which an
 * earlier build registered. Safe to call on every load.
 */
export async function removeLegacyDownloadWorker(): Promise<void> {
  if (typeof navigator === 'undefined' || !('serviceWorker' in navigator)) return;
  const script = new URL(WORKER_FILE, document.baseURI).href;
  const scope = new URL(DOWNLOAD_PATH, document.baseURI).href;
  const registrations = await navigator.serviceWorker.getRegistrations();
  await Promise.all(
    registrations
      .filter(
        (registration) =>
          registration.scope !== scope &&
          [registration.active, registration.waiting, registration.installing].some(
            (worker) => worker?.scriptURL === script,
          ),
      )
      .map((registration) => registration.unregister()),
  );
}

export class StreamStartError extends Error {}

/**
 * Save `source` as `filename` through the download worker. Rejects with a
 * {@link StreamStartError} if the download could not begin (the caller may
 * then fall back), or with the source's error if it fails partway.
 */
export async function streamDownload({
  filename,
  source,
  signal,
}: {
  filename: string;
  source: AsyncIterable<Uint8Array>;
  signal?: AbortSignal;
}): Promise<void> {
  let registration: ServiceWorkerRegistration;
  try {
    registration = await downloadWorker();
  } catch (error) {
    throw new StreamStartError(error instanceof Error ? error.message : String(error));
  }
  const worker = registration.active;
  if (!worker) throw new StreamStartError('The download helper is not active.');

  const channel = new MessageChannel();
  const server = servePort(channel.port1 as unknown as PortLike, source, { signal });
  const id = crypto.randomUUID();
  worker.postMessage({ type: 'register', id, filename }, [channel.port2]);

  const iframe = document.createElement('iframe');
  iframe.hidden = true;
  const ping = setInterval(
    () => registration.active?.postMessage({ type: 'ping' }),
    PING_MS,
  );
  try {
    try {
      await Promise.race([
        server.registered,
        server.done,
        timeout(START_TIMEOUT_MS, 'The download helper did not answer.'),
      ]);
      // The browser may have stopped and restarted the worker since, losing
      // its in-memory list, so register again on a fresh channel right before
      // the browser asks for the URL.
      const again = new MessageChannel();
      const reRegistered = server.adopt(again.port1 as unknown as PortLike);
      (registration.active ?? worker).postMessage({ type: 'register', id, filename }, [
        again.port2,
      ]);
      await Promise.race([
        reRegistered,
        server.done,
        timeout(START_TIMEOUT_MS, 'The download helper did not answer.'),
      ]);
      iframe.src = downloadUrl(registration.scope, id, filename);
      document.body.appendChild(iframe);
      await Promise.race([
        server.started,
        server.done,
        timeout(START_TIMEOUT_MS, 'The browser did not start the download.'),
      ]);
    } catch (error) {
      // Tells the worker to drop the download and stops the source.
      server.stop(error);
      if (error instanceof DOMException && error.name === 'AbortError') throw error;
      throw new StreamStartError(
        error instanceof Error ? error.message : String(error),
      );
    }
    await server.done;
  } finally {
    clearInterval(ping);
    // Give the browser a moment to hand the response to its download manager.
    setTimeout(() => iframe.remove(), 1_000);
  }
}
