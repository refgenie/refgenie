/**
 * Shared FASTA digest module. Source of truth: refget, frontend/src/lib/fastaDigest/.
 * Other apps (refgenie) copy this folder verbatim; change it in refget, then copy.
 *
 * Compute a FASTA file's GA4GH sequence collection digest in the browser.
 *
 * The file is read in chunks and fed to the `@databio/gtars` wasm streaming
 * hasher (`fastaHasherNew` / `Update` / `Finish`), so memory stays small
 * whatever the file size. Gzip is decompressed here, not in wasm: the wasm
 * decoder in @databio/gtars 0.9.2 reads only the first gzip member, and the
 * browser's `DecompressionStream` rejects a second one, so a bgzip (BGZF)
 * file, which is a run of small gzip members, is split into its members by
 * their recorded sizes and each member is decompressed on its own.
 *
 * No React or DOM here, so it runs under Node tests. The browser runs it in a
 * Web Worker (`fastaDigest.worker.ts`), started by `fastaDigestClient.ts`.
 */

type Gtars = typeof import('@databio/gtars');
export type FastaHasherWasm = Pick<
  Gtars,
  'fastaHasherNew' | 'fastaHasherUpdate' | 'fastaHasherFinish' | 'fastaHasherFree'
>;

/** One sequence of the file, as the wasm hasher reports it. */
export interface FastaDigestSequence {
  name: string;
  length: number;
  sha512t24u: string;
  md5: string;
  alphabet: string;
  description: string | null;
}

/** What the wasm hasher returns for one file. */
export interface FastaDigestResult {
  digest: string;
  n_sequences: number;
  names_digest: string;
  sequences_digest: string;
  lengths_digest: string;
  sequences: FastaDigestSequence[];
}

export type Compression = 'none' | 'gzip' | 'bgzf';

/** Bytes always read before telling plain text, gzip and BGZF apart. */
export const SNIFF_BYTES = 64;

/**
 * How many bytes `sniffCompression` needs to see to decide about `head`: at
 * least `SNIFF_BYTES`, and for a gzip header with an extra field, the whole
 * extra field (12 + XLEN bytes), which may be longer.
 */
export function sniffBytesNeeded(head: Uint8Array): number {
  if (head.length < 12 || head[0] !== 0x1f || head[1] !== 0x8b) return SNIFF_BYTES;
  if ((head[3] & 0x04) === 0) return SNIFF_BYTES;
  const xlen = head[10] | (head[11] << 8);
  return Math.max(SNIFF_BYTES, 12 + xlen);
}

/**
 * Plain text, ordinary gzip, or BGZF: a gzip member whose extra field holds a
 * `BC` subfield (the member's own size).
 */
export function sniffCompression(head: Uint8Array): Compression {
  if (head.length < 2 || head[0] !== 0x1f || head[1] !== 0x8b) return 'none';
  return bgzfBlockSize(head, 0) === null ? 'gzip' : 'bgzf';
}

/**
 * Total size of the BGZF member at `offset` (header to trailer), or `null`
 * when the bytes there are not a BGZF header or are too short to tell.
 */
export function bgzfBlockSize(buf: Uint8Array, offset: number): number | null {
  if (buf.length - offset < 12) return null;
  if (buf[offset] !== 0x1f || buf[offset + 1] !== 0x8b || buf[offset + 2] !== 8)
    return null;
  if ((buf[offset + 3] & 0x04) === 0) return null; // no FEXTRA
  const xlen = buf[offset + 10] | (buf[offset + 11] << 8);
  const end = offset + 12 + xlen;
  if (buf.length < end) return null;
  for (let p = offset + 12; p + 4 <= end;) {
    const slen = buf[p + 2] | (buf[p + 3] << 8);
    if (buf[p] === 0x42 && buf[p + 1] === 0x43 && slen === 2 && p + 6 <= end) {
      return (buf[p + 4] | (buf[p + 5] << 8)) + 1;
    }
    p += 4 + slen;
  }
  return null;
}

/** Concatenate `a` and `b`, reusing `b` when `a` is empty. */
function join(a: Uint8Array, b: Uint8Array): Uint8Array {
  if (a.length === 0) return b;
  const out = new Uint8Array(a.length + b.length);
  out.set(a);
  out.set(b, a.length);
  return out;
}

export const BGZF_CUT_OFF =
  'This bgzip file looks cut off: it has no end-of-file marker.';

export const NOT_BGZF =
  'This gzip file is made of several parts that are not bgzip blocks. Decompress it, or compress it again with bgzip or gzip, and try again.';

/**
 * Split a BGZF byte stream into its whole gzip members. A bgzip file ends with
 * an empty member (the EOF marker); a file whose last member is not empty was
 * cut off at a block boundary, so this throws rather than yield a short file.
 */
export async function* bgzfMembers(
  chunks: AsyncIterable<Uint8Array>,
): AsyncGenerator<Uint8Array> {
  let buf: Uint8Array = new Uint8Array(0);
  let pos = 0;
  let lastWasEof = false;
  for await (const chunk of chunks) {
    buf = join(buf.subarray(pos), chunk);
    pos = 0;
    for (;;) {
      if (buf.length - pos < 12) break;
      const size = bgzfBlockSize(buf, pos);
      if (size === null) {
        // A whole header is here but it is not BGZF: the file mixes formats.
        const xlen = buf[pos + 10] | (buf[pos + 11] << 8);
        if (buf.length - pos >= 12 + xlen) throw new Error(NOT_BGZF);
        break;
      }
      if (buf.length - pos < size) break;
      const member = buf.subarray(pos, pos + size);
      // ISIZE, the member's uncompressed size, is its last 4 bytes.
      lastWasEof =
        (member[size - 4] | member[size - 3] | member[size - 2] | member[size - 1]) === 0;
      yield member;
      pos += size;
    }
  }
  if (buf.length - pos > 0)
    throw new Error('The file ends partway through a gzip block.');
  if (!lastWasEof) throw new Error(BGZF_CUT_OFF);
}

export type Gunzip = (member: Uint8Array) => Promise<Uint8Array>;

/** Decompress one whole gzip member with the browser's `DecompressionStream`. */
export const gunzipWeb: Gunzip = async (member) => {
  const stream = new Blob([member as BlobPart])
    .stream()
    .pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
};

/**
 * Decompress members `concurrency` at a time and yield their bytes in order.
 * BGZF members hold at most 64 KiB, so a batch is small.
 */
export async function* inflateMembers(
  members: AsyncIterable<Uint8Array>,
  gunzip: Gunzip = gunzipWeb,
  concurrency = 16,
): AsyncGenerator<Uint8Array> {
  let batch: Promise<Uint8Array>[] = [];
  for await (const member of members) {
    const out = gunzip(member);
    // If the consumer stops early, the rest of the batch is never awaited.
    out.catch(() => {});
    batch.push(out);
    if (batch.length >= concurrency) {
      for (const out of await Promise.all(batch)) yield out;
      batch = [];
    }
  }
  for (const out of await Promise.all(batch)) yield out;
}

/** Stream-decompress an ordinary (single-member) gzip file. */
export async function* gunzipStream(
  chunks: AsyncIterable<Uint8Array>,
): AsyncGenerator<Uint8Array> {
  const ds = new DecompressionStream('gzip');
  const writer = ds.writable.getWriter();
  const pump = (async () => {
    try {
      for await (const chunk of chunks) {
        await writer.write(chunk as BufferSource);
      }
      await writer.close();
    } catch (err) {
      await writer.abort(err).catch(() => {});
    }
  })();
  const reader = ds.readable.getReader();
  let done = false;
  try {
    for (;;) {
      const next = await reader.read();
      if (next.done) break;
      yield next.value;
    }
    done = true;
  } catch (err) {
    const message = err instanceof Error ? err.message : String(err);
    if (/junk|trailing/i.test(message)) {
      throw Object.assign(new Error(NOT_BGZF), { cause: err });
    }
    throw err;
  } finally {
    // Stopped early: cancel so the pump's pending write settles.
    if (!done) await reader.cancel().catch(() => {});
    await pump;
  }
}

/** Read a Blob as a stream of chunks. */
export async function* blobChunks(blob: Blob): AsyncGenerator<Uint8Array> {
  const reader = blob.stream().getReader();
  try {
    for (;;) {
      const { done, value } = await reader.read();
      if (done) return;
      yield value;
    }
  } finally {
    reader.releaseLock();
  }
}

/** Prepend `head` to the rest of a chunk stream. */
async function* prepend(
  head: Uint8Array,
  rest: AsyncIterator<Uint8Array>,
): AsyncGenerator<Uint8Array> {
  let done = false;
  try {
    if (head.length > 0) yield head;
    for (;;) {
      const next = await rest.next();
      if (next.done) {
        done = true;
        return;
      }
      yield next.value;
    }
  } finally {
    // Stopped early: close the upstream so it releases its reader.
    if (!done) await rest.return?.();
  }
}

export interface DigestOptions {
  /** Called as file (compressed) bytes are read: `(bytesRead)`. */
  onProgress?: (bytesRead: number) => void;
  gunzip?: Gunzip;
}

/**
 * Digest a FASTA (plain, gzip or bgzip) given as a chunk stream. Throws on a
 * file the hasher cannot parse.
 */
export async function digestFastaChunks(
  chunks: AsyncIterable<Uint8Array>,
  wasm: FastaHasherWasm,
  { onProgress, gunzip }: DigestOptions = {},
): Promise<FastaDigestResult> {
  let read = 0;
  const counted = (async function* () {
    for await (const chunk of chunks) {
      read += chunk.length;
      onProgress?.(read);
      yield chunk;
    }
  })();

  // Sniff the format from the first bytes.
  const it = counted[Symbol.asyncIterator]();
  let head: Uint8Array = new Uint8Array(0);
  while (head.length < sniffBytesNeeded(head)) {
    const next = await it.next();
    if (next.done) break;
    head = join(head, next.value);
  }
  const compression = sniffCompression(head);
  const raw = prepend(head, it);
  const plain =
    compression === 'bgzf'
      ? inflateMembers(bgzfMembers(raw), gunzip)
      : compression === 'gzip'
        ? gunzipStream(raw)
        : raw;

  const handle = wasm.fastaHasherNew();
  let finished = false;
  try {
    let sawData = false;
    for await (const chunk of plain) {
      if (chunk.length === 0) continue;
      if (!sawData && chunk[0] === 0x1f && chunk[1] === 0x8b) {
        throw new Error(
          'This file is compressed twice. Decompress it once and try again.',
        );
      }
      sawData = true;
      wasm.fastaHasherUpdate(handle, chunk);
    }
    if (!sawData) throw new Error('The file is empty.');
    finished = true;
    const result = wasm.fastaHasherFinish(handle) as FastaDigestResult;
    if (!result || result.n_sequences === 0) {
      throw new Error('No sequences found. Is this a FASTA file?');
    }
    return result;
  } finally {
    if (!finished) wasm.fastaHasherFree(handle);
  }
}

/** File names this module accepts (for an `<input accept>`). */
export const FASTA_ACCEPT =
  '.fa,.fasta,.fna,.fas,.fa.gz,.fasta.gz,.fna.gz,.fas.gz,.gz,.bgz';

/** True when the name looks like a FASTA file (optionally gzipped). */
export function looksLikeFasta(name: string): boolean {
  return /\.(fa|fasta|fna|fas)(\.b?gz)?$/i.test(name);
}
