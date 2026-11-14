/**
 * Shared FASTA digest module. Source of truth: refget, frontend/src/lib/fastaDigest/.
 * Other apps (refgenie) copy this folder verbatim; change it in refget, then copy.
 *
 * The page's side of the FASTA digest worker: start a digest, hear status and
 * progress, cancel. Tests swap the worker factory for a fake.
 */

import type { FastaDigestResult } from './fastaDigest.ts';

export interface DigestWorkerRequest {
  file: File;
}

/** Timing counters for one digest. `totalBytes` is the file's size on disk. */
export interface DigestStats {
  chunks: number;
  totalBytes: number;
  elapsedMs: number;
  avgChunkSize: number;
}

export type DigestWorkerReply =
  | { type: 'status'; message: string }
  | { type: 'progress'; bytesRead: number; totalBytes: number }
  | { type: 'result'; result: FastaDigestResult; stats: DigestStats }
  | { type: 'error'; message: string };

/** What a finished digest resolves to. */
export interface DigestOutcome {
  result: FastaDigestResult;
  stats: DigestStats;
}

export interface DigestCallbacks {
  /** File bytes read so far (compressed bytes for a gzipped file). */
  onProgress?: (bytesRead: number, totalBytes: number) => void;
  /** A short line saying what the worker is doing. */
  onStatus?: (message: string) => void;
}

/** The slice of `Worker` used here. */
export interface WorkerLike {
  postMessage(message: DigestWorkerRequest): void;
  onmessage: ((event: MessageEvent<DigestWorkerReply>) => void) | null;
  onerror: ((event: ErrorEvent) => void) | null;
  terminate(): void;
}

export type WorkerFactory = () => WorkerLike;

export const createDigestWorker: WorkerFactory = () =>
  new Worker(new URL('./fastaDigest.worker.ts', import.meta.url), {
    type: 'module',
  }) as unknown as WorkerLike;

export interface DigestJob {
  done: Promise<DigestOutcome>;
  /** Stop at once; `done` rejects with an AbortError. */
  cancel: () => void;
}

/** True for the rejection `cancel` causes. */
export function isAbortError(err: unknown): boolean {
  return err instanceof DOMException && err.name === 'AbortError';
}

/** Digest `file` in a worker. */
export function digestFile(
  file: File,
  { onProgress, onStatus }: DigestCallbacks = {},
  factory: WorkerFactory = createDigestWorker,
): DigestJob {
  const worker = factory();
  let settle: { reject: (err: unknown) => void } | null = null;
  const done = new Promise<DigestOutcome>((resolve, reject) => {
    settle = { reject };
    worker.onmessage = (event) => {
      const reply = event.data;
      if (reply.type === 'status') onStatus?.(reply.message);
      else if (reply.type === 'progress') onProgress?.(reply.bytesRead, reply.totalBytes);
      else if (reply.type === 'result') {
        worker.terminate();
        resolve({ result: reply.result, stats: reply.stats });
      } else {
        worker.terminate();
        reject(new Error(reply.message));
      }
    };
    worker.onerror = (event) => {
      event.preventDefault?.();
      worker.terminate();
      reject(new Error(event.message || 'The digest worker failed to start.'));
    };
  });
  worker.postMessage({ file });
  return {
    done,
    cancel: () => {
      worker.terminate();
      settle?.reject(new DOMException('Cancelled', 'AbortError'));
    },
  };
}
