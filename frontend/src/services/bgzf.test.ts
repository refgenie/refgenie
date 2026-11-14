// @vitest-environment node
import { gunzipSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { BGZF_BLOCK_INPUT, BGZF_EOF, bgzf, bgzfBlock, crc32 } from './bgzf';

async function collect(source: AsyncIterable<Uint8Array>): Promise<Uint8Array> {
  const parts: Uint8Array[] = [];
  for await (const part of source) parts.push(part);
  return Buffer.concat(parts);
}

async function* chunked(data: Uint8Array, sizes: number[]) {
  let offset = 0;
  let i = 0;
  while (offset < data.length) {
    const size = sizes[i++ % sizes.length];
    yield data.slice(offset, offset + size);
    offset += size;
  }
}

/** Split a BGZF file into members using each header's BSIZE. */
function members(file: Uint8Array): Uint8Array[] {
  const out: Uint8Array[] = [];
  const view = new DataView(file.buffer, file.byteOffset, file.byteLength);
  let offset = 0;
  while (offset < file.length) {
    expect([
      file[offset],
      file[offset + 1],
      file[offset + 2],
      file[offset + 3],
    ]).toEqual([0x1f, 0x8b, 8, 4]);
    expect(view.getUint16(offset + 10, true)).toBe(6);
    expect(String.fromCharCode(file[offset + 12], file[offset + 13])).toBe('BC');
    const size = view.getUint16(offset + 16, true) + 1;
    out.push(file.subarray(offset, offset + size));
    offset += size;
  }
  expect(offset).toBe(file.length);
  return out;
}

function fakeFasta(bytes: number): Uint8Array {
  const bases = 'ACGTTGCANNACGT';
  const out = new Uint8Array(bytes);
  let seed = 7;
  for (let i = 0; i < bytes; i += 1) {
    seed = (seed * 1103515245 + 12345) & 0x7fffffff;
    out[i] = i % 81 === 80 ? 10 : bases.charCodeAt(seed % bases.length);
  }
  return out;
}

describe('crc32', () => {
  it('matches the standard check value', () => {
    expect(crc32(new TextEncoder().encode('123456789'))).toBe(0xcbf43926);
  });
});

describe('bgzf', () => {
  it('gunzips to the exact input, with correct BSIZE in every member', async () => {
    const source = fakeFasta(BGZF_BLOCK_INPUT * 3 + 1234);
    const file = await collect(bgzf(chunked(source, [1000, 70_000, 17])));
    expect(Buffer.from(gunzipSync(file)).equals(Buffer.from(source))).toBe(true);
    const blocks = members(file);
    expect(blocks).toHaveLength(5); // 4 data blocks + EOF
    expect([...blocks[4]]).toEqual([...BGZF_EOF]);
    for (const block of blocks.slice(0, 4)) {
      expect(block.length).toBeLessThanOrEqual(65536);
      // Each member is a valid gzip file on its own.
      expect(gunzipSync(block).length).toBeLessThanOrEqual(BGZF_BLOCK_INPUT);
    }
  });

  it('writes only the EOF block for empty input', async () => {
    const file = await collect(bgzf(chunked(new Uint8Array(0), [1])));
    expect([...file]).toEqual([...BGZF_EOF]);
    expect(gunzipSync(file).length).toBe(0);
  });

  it('stores a block that would not fit compressed', async () => {
    const data = new Uint8Array(BGZF_BLOCK_INPUT).fill(65);
    // A deflater that makes things bigger forces the stored-block path.
    const block = await bgzfBlock(data, async () => new Uint8Array(70_000));
    expect(block.length).toBe(18 + 5 + data.length + 8);
    expect(Buffer.from(gunzipSync(block)).equals(Buffer.from(data))).toBe(true);
  });

  it('writes no EOF block when the source fails', async () => {
    async function* failing() {
      yield fakeFasta(BGZF_BLOCK_INPUT * 2);
      throw new Error('boom');
    }
    const seen: Uint8Array[] = [];
    await expect(
      (async () => {
        for await (const part of bgzf(failing())) seen.push(part);
      })(),
    ).rejects.toThrow('boom');
    expect(seen.some((part) => Buffer.from(part).equals(Buffer.from(BGZF_EOF)))).toBe(
      false,
    );
  });
});
