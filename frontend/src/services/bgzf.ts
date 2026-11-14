/**
 * BGZF (blocked gzip) compression in the browser, as written by `bgzip`.
 *
 * A BGZF file is a series of gzip members, each holding at most 64 KiB of
 * input and recording its own compressed size in a `BC` extra field, then a
 * fixed 28-byte empty member that marks the end. Plain gzip tools read it as
 * one gzip file; `samtools faidx` can index it and seek within it.
 *
 * Each block is deflated with the browser's `CompressionStream('deflate-raw')`;
 * the header, CRC32 and size trailer are written here. No DOM, so it runs
 * under vitest in Node.
 */

/** Input bytes per block: bgzip's own choice, so the compressed block fits in 64 KiB. */
export const BGZF_BLOCK_INPUT = 0xff00;
/** Largest whole member, header and trailer included. */
const MAX_BLOCK_SIZE = 0x10000;
const HEADER_SIZE = 18;
const TRAILER_SIZE = 8;

/** The empty member that ends every BGZF file. */
export const BGZF_EOF = Uint8Array.from([
  0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
  0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
]);

const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n += 1) {
    let c = n;
    for (let k = 0; k < 8; k += 1) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

/** The gzip CRC-32 of `data`. */
export function crc32(data: Uint8Array): number {
  let crc = 0xffffffff;
  for (let i = 0; i < data.length; i += 1) {
    crc = CRC_TABLE[(crc ^ data[i]) & 0xff] ^ (crc >>> 8);
  }
  return (crc ^ 0xffffffff) >>> 0;
}

/** Raw-deflate one buffer (no zlib or gzip wrapper). */
export type DeflateRaw = (data: Uint8Array) => Promise<Uint8Array>;

/** {@link DeflateRaw} with the Compression Streams API. */
export const deflateRawWeb: DeflateRaw = async (data) => {
  const stream = new CompressionStream('deflate-raw');
  const writer = stream.writable.getWriter();
  // Not awaited: the write only finishes as the readable side is drained.
  writer.write(data as Uint8Array<ArrayBuffer>).catch(() => {});
  writer.close().catch(() => {});
  return new Uint8Array(await new Response(stream.readable).arrayBuffer());
};

/** `data` as one uncompressed ("stored") deflate block. */
function storedDeflate(data: Uint8Array): Uint8Array {
  const out = new Uint8Array(data.length + 5);
  out[0] = 0x01; // final block, stored
  out[1] = data.length & 0xff;
  out[2] = data.length >>> 8;
  out[3] = ~data.length & 0xff;
  out[4] = (~data.length >>> 8) & 0xff;
  out.set(data, 5);
  return out;
}

/** One BGZF member holding `data` (at most {@link BGZF_BLOCK_INPUT} bytes). */
export async function bgzfBlock(data: Uint8Array, deflate: DeflateRaw = deflateRawWeb) {
  if (data.length > BGZF_BLOCK_INPUT) throw new Error('BGZF block input too large');
  let body = await deflate(data);
  // Data that does not compress could overflow the block; store it instead.
  if (body.length + HEADER_SIZE + TRAILER_SIZE > MAX_BLOCK_SIZE)
    body = storedDeflate(data);
  const size = HEADER_SIZE + body.length + TRAILER_SIZE;
  const out = new Uint8Array(size);
  const view = new DataView(out.buffer);
  // gzip header: magic, deflate, FEXTRA, no mtime, no extra flags, unknown OS.
  out.set([0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff], 0);
  view.setUint16(10, 6, true); // XLEN
  out[12] = 0x42; // 'B'
  out[13] = 0x43; // 'C'
  view.setUint16(14, 2, true); // SLEN
  view.setUint16(16, size - 1, true); // BSIZE: member size minus one
  out.set(body, HEADER_SIZE);
  view.setUint32(HEADER_SIZE + body.length, crc32(data), true);
  view.setUint32(HEADER_SIZE + body.length + 4, data.length, true);
  return out;
}

export interface BgzfOptions {
  deflate?: DeflateRaw;
  /** Blocks compressed at once. */
  concurrency?: number;
}

/**
 * Compress a byte stream to BGZF, one member per 64 KiB of input, in order.
 * The end-of-file block is written only once the source finishes without
 * error, so a failed run never looks like a complete file.
 */
export async function* bgzf(
  source: AsyncIterable<Uint8Array>,
  { deflate = deflateRawWeb, concurrency = 4 }: BgzfOptions = {},
): AsyncGenerator<Uint8Array> {
  const queue: Promise<Uint8Array>[] = [];
  let buffer = new Uint8Array(BGZF_BLOCK_INPUT);
  let filled = 0;

  const enqueue = (data: Uint8Array) => {
    const pending = bgzfBlock(data, deflate);
    // Awaited in order below; this only keeps an early failure from going unhandled.
    pending.catch(() => {});
    queue.push(pending);
  };

  for await (const chunk of source) {
    let offset = 0;
    while (offset < chunk.length) {
      const take = Math.min(BGZF_BLOCK_INPUT - filled, chunk.length - offset);
      buffer.set(chunk.subarray(offset, offset + take), filled);
      filled += take;
      offset += take;
      if (filled === BGZF_BLOCK_INPUT) {
        enqueue(buffer);
        buffer = new Uint8Array(BGZF_BLOCK_INPUT);
        filled = 0;
        if (queue.length >= concurrency)
          yield await (queue.shift() as Promise<Uint8Array>);
      }
    }
  }
  if (filled > 0) enqueue(buffer.subarray(0, filled));
  while (queue.length > 0) yield await (queue.shift() as Promise<Uint8Array>);
  yield BGZF_EOF;
}
