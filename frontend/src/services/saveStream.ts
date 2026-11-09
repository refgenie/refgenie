/**
 * Save a stream of bytes to a file, without holding the whole file in memory
 * where the browser allows it.
 *
 * - Chrome and Edge: the File System Access API writes straight to the file
 *   the user picks.
 * - Firefox and Safari: a service worker streams the bytes to the browser's
 *   download manager (`streamDownload.ts`).
 * - Otherwise (no service worker, e.g. a plain http:// page that is not
 *   localhost): the bytes are collected into a Blob, so only small files.
 */

import { StreamStartError, canStreamDownload, streamDownload } from './streamDownload';

/** Largest file the Blob fallback will build in memory. */
export const BLOB_FALLBACK_MAX = 500_000_000;

export class FileTooLargeError extends Error {
  constructor() {
    super(
      'This file is too large to save from this page in this browser. ' +
        'Choose fewer sequences, or use the command line.',
    );
    this.name = 'FileTooLargeError';
  }
}

export type SaveMethod = 'file-picker' | 'stream' | 'blob';

/**
 * How to save a file of about `bytes`: the File System Access API if present,
 * else a service-worker stream, else a Blob if the file is small. `null` if
 * none can do it.
 */
export function chooseSaveMethod({
  hasFilePicker,
  canStream,
  bytes,
}: {
  hasFilePicker: boolean;
  canStream: boolean;
  bytes: number;
}): SaveMethod | null {
  if (hasFilePicker) return 'file-picker';
  if (canStream) return 'stream';
  return bytes <= BLOB_FALLBACK_MAX ? 'blob' : null;
}

interface WritableFileStream {
  write(data: Uint8Array): Promise<void>;
  close(): Promise<void>;
  abort(): Promise<void>;
}

interface SaveFileHandle {
  createWritable(): Promise<WritableFileStream>;
}

type ShowSaveFilePicker = (options: {
  suggestedName: string;
  types: Array<{ description: string; accept: Record<string, string[]> }>;
}) => Promise<SaveFileHandle>;

export interface SaveStreamOptions {
  suggestedName: string;
  /** Size estimate in bytes, to rule out a Blob that would be too large. */
  estimatedBytes: number;
  /** The picker's file type, e.g. `{ description: 'FASTA', accept: { 'application/gzip': ['.gz'] } }`. */
  fileType: { description: string; accept: Record<string, string[]> };
  chunks: () => AsyncGenerator<Uint8Array>;
  signal?: AbortSignal;
}

/**
 * Write the chunks to a file. Call it straight from a click handler: the save
 * picker must open while the user gesture is still active.
 */
export async function saveStream({
  suggestedName,
  estimatedBytes,
  fileType,
  chunks,
  signal,
}: SaveStreamOptions): Promise<void> {
  const win = window as unknown as { showSaveFilePicker?: ShowSaveFilePicker };
  const method = chooseSaveMethod({
    hasFilePicker: !!win.showSaveFilePicker,
    canStream: canStreamDownload(),
    bytes: estimatedBytes,
  });

  if (method === 'file-picker') {
    // Called on `window`: an unbound call throws "Illegal invocation".
    const handle = await win.showSaveFilePicker!({ suggestedName, types: [fileType] });
    const writable = await handle.createWritable();
    try {
      for await (const chunk of chunks()) await writable.write(chunk);
      await writable.close();
    } catch (error) {
      await writable.abort().catch(() => {});
      throw error;
    }
    return;
  }

  if (method === 'stream') {
    try {
      await streamDownload({ filename: suggestedName, source: chunks(), signal });
      return;
    } catch (error) {
      // Only a download that never began can fall back.
      if (!(error instanceof StreamStartError)) throw error;
      if (estimatedBytes > BLOB_FALLBACK_MAX) {
        throw new Error(`Could not start the download: ${error.message}`);
      }
    }
  }

  if (method === null) throw new FileTooLargeError();
  await saveBlob(suggestedName, chunks());
}

async function saveBlob(
  filename: string,
  source: AsyncIterable<Uint8Array>,
): Promise<void> {
  const parts: BlobPart[] = [];
  for await (const chunk of source) parts.push(chunk as Uint8Array<ArrayBuffer>);
  const url = URL.createObjectURL(
    new Blob(parts, { type: 'application/octet-stream' }),
  );
  try {
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
  } finally {
    // Revoke after the click has started the download.
    setTimeout(() => URL.revokeObjectURL(url), 0);
  }
}
