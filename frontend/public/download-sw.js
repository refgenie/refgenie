/*
 * Streaming download service worker.
 *
 * Browsers without the File System Access API (Firefox, Safari) can still save
 * a large file without holding it in memory: the page registers a download
 * here over a MessageChannel, then opens a unique URL under this worker's
 * scope in a hidden iframe. The scope is the download path itself
 * (`<base>/rg-download/`), so the worker controls no other page. This worker answers that URL with a streamed
 * response marked as an attachment, and pulls each chunk from the page over
 * the channel only when the browser asks for more, so backpressure holds.
 *
 * Messages, page -> worker (on the worker): {type: 'register', id, filename}
 * with a MessagePort (sent again, with a new port, just before the iframe
 * opens; the newest one wins), and {type: 'ping'} to keep the worker alive.
 * A registration the browser never fetches is dropped when the page sends
 * 'error' on its port, or after EXPIRE_MS.
 * On the port, worker -> page: 'registered', 'started', 'pull', 'cancel';
 * page -> worker: {type: 'chunk', chunk}, {type: 'end'}, {type: 'error', message}.
 *
 * Keep in step with src/services/streamDownload.ts.
 */

/* global self */

const PREFIX = 'rg-download/';
const EXPIRE_MS = 60_000;
const pending = new Map();

/** The path this worker answers: its scope, or `rg-download/` under an older, wider one. */
function downloadRoot() {
  const scope = new URL(self.registration.scope).pathname;
  return scope.endsWith('/' + PREFIX) ? scope : scope + PREFIX;
}

/** Drop a registration that was never fetched. */
function forget(id, port) {
  const entry = pending.get(id);
  if (entry && entry.port === port) {
    clearTimeout(entry.timer);
    pending.delete(id);
  }
  port.onmessage = null;
  port.close();
}

function downloadStream(port) {
  let controller = null;
  let wake = null;
  let finished = false;

  const stream = new ReadableStream(
    {
      start(c) {
        controller = c;
      },
      pull() {
        return new Promise((resolve) => {
          wake = resolve;
          port.postMessage({ type: 'pull' });
        });
      },
      cancel() {
        // The user cancelled in the browser's download list.
        finished = true;
        port.postMessage({ type: 'cancel' });
        port.close();
      },
    },
    { highWaterMark: 2 },
  );

  port.onmessage = (event) => {
    if (finished) return;
    const msg = event.data || {};
    if (msg.type === 'chunk') {
      controller.enqueue(msg.chunk);
    } else if (msg.type === 'end') {
      finished = true;
      controller.close();
      port.close();
    } else if (msg.type === 'error') {
      // Error the stream, so the browser marks the download failed instead of
      // keeping a truncated file that looks complete.
      finished = true;
      controller.error(new Error(msg.message || 'Download failed'));
      port.close();
    }
    const resolve = wake;
    wake = null;
    if (resolve) resolve();
  };

  return stream;
}

self.addEventListener('install', () => {
  self.skipWaiting();
});

self.addEventListener('activate', (event) => {
  event.waitUntil(self.clients.claim());
});

self.addEventListener('message', (event) => {
  const msg = event.data || {};
  if (msg.type === 'register' && event.ports && event.ports[0]) {
    const port = event.ports[0];
    const id = msg.id;
    const old = pending.get(id);
    if (old) forget(id, old.port);
    // Until the browser fetches the URL, the page only ever says it gave up.
    port.onmessage = (e) => {
      const m = e.data || {};
      if (m.type === 'error' || m.type === 'end') forget(id, port);
    };
    const timer = setTimeout(() => forget(id, port), EXPIRE_MS);
    pending.set(id, { port, filename: msg.filename, timer });
    port.postMessage({ type: 'registered' });
  }
  // 'ping' needs no answer: receiving it keeps the worker alive.
});

self.addEventListener('fetch', (event) => {
  const url = new URL(event.request.url);
  const root = downloadRoot();
  if (!url.pathname.startsWith(root)) return;
  const id = url.pathname.slice(root.length).split('/')[0];
  const entry = pending.get(id);
  if (!entry) {
    event.respondWith(new Response('This download has expired.', { status: 404 }));
    return;
  }
  clearTimeout(entry.timer);
  pending.delete(id);
  const stream = downloadStream(entry.port);
  entry.port.postMessage({ type: 'started' });
  const filename = entry.filename;
  event.respondWith(
    new Response(stream, {
      headers: {
        'Content-Type': 'application/octet-stream',
        'Content-Disposition':
          `attachment; filename="${filename.replace(/["\\\r\n]/g, '_')}"; ` +
          `filename*=UTF-8''${encodeURIComponent(filename)}`,
        'X-Content-Type-Options': 'nosniff',
      },
    }),
  );
});
