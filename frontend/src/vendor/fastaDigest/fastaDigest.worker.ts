/**
 * Shared FASTA digest module. Source of truth: refget, frontend/src/lib/fastaDigest/.
 * Other apps (refgenie) copy this folder verbatim; change it in refget, then copy.
 *
 * Web Worker: digest one FASTA file off the main thread, so the page stays
 * responsive on a 3 GB genome. The page cancels by terminating the worker.
 * Keep the message shapes in step with `fastaDigestClient.ts`.
 */

import { blobChunks, digestFastaChunks } from './fastaDigest.ts';
import type { DigestWorkerReply, DigestWorkerRequest } from './fastaDigestClient.ts';

// The DOM and WebWorker type libraries cannot both be loaded into one
// program, so name the slice of the worker global used here.
interface WorkerScope {
  onmessage: ((event: MessageEvent<DigestWorkerRequest>) => void) | null;
  postMessage(message: DigestWorkerReply): void;
}
const ctx = self as unknown as WorkerScope;

const PROGRESS_MS = 150;

ctx.onmessage = async (event) => {
  const { file } = event.data;
  const startTime = Date.now();
  let chunks = 0;
  async function* counted(): AsyncGenerator<Uint8Array> {
    for await (const chunk of blobChunks(file)) {
      chunks++;
      yield chunk;
    }
  }
  try {
    ctx.postMessage({ type: 'status', message: 'Loading WASM module...' });
    const mod = await import('@databio/gtars');
    await mod.default();
    ctx.postMessage({ type: 'status', message: 'Processing file...' });
    let last = 0;
    const result = await digestFastaChunks(counted(), mod, {
      onProgress: (bytesRead) => {
        const now = Date.now();
        if (now - last < PROGRESS_MS) return;
        last = now;
        ctx.postMessage({ type: 'progress', bytesRead, totalBytes: file.size });
      },
    });
    ctx.postMessage({ type: 'progress', bytesRead: file.size, totalBytes: file.size });
    const elapsedMs = Date.now() - startTime;
    ctx.postMessage({
      type: 'result',
      result,
      stats: {
        chunks,
        totalBytes: file.size,
        elapsedMs,
        avgChunkSize: chunks > 0 ? Math.round(file.size / chunks) : 0,
      },
    });
  } catch (err) {
    ctx.postMessage({
      type: 'error',
      message: err instanceof Error ? err.message : String(err),
    });
  }
};
