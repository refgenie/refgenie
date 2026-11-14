// @vitest-environment node
/// <reference types="node" />
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { gzipSync } from 'node:zlib';
import { beforeAll, describe, expect, it } from 'vitest';
import { BGZF_EOF, bgzfBlock } from './bgzf';
import {
  bgzfBlockSize,
  bgzfMembers,
  digestFastaChunks,
  looksLikeFasta,
  sniffCompression,
  type FastaHasherWasm,
} from '../vendor/fastaDigest/fastaDigest.ts';

const FIXTURE = join(
  dirname(fileURLToPath(import.meta.url)),
  '../test/fixtures/fasta-store',
);
const SOURCE = new Uint8Array(readFileSync(join(FIXTURE, 'source.fa')));
// Computed by the refget CLI when the fixture store was built.
const DIGEST = readFileSync(join(FIXTURE, 'digest.txt'), 'utf8').trim();

let wasm: FastaHasherWasm;

beforeAll(async () => {
  const dir = dirname(createRequire(import.meta.url).resolve('@databio/gtars'));
  const mod = await import('@databio/gtars');
  await mod.default({ module_or_path: readFileSync(join(dir, 'gtars_js_bg.wasm')) });
  wasm = mod;
});

/** `bytes` as a stream of `size`-byte chunks. */
async function* chunked(bytes: Uint8Array, size: number): AsyncGenerator<Uint8Array> {
  for (let i = 0; i < bytes.length; i += size) yield bytes.slice(i, i + size);
}

/** BGZF with small blocks, so a short file still spans several members. */
async function bgzipped(bytes: Uint8Array, blockInput = 300): Promise<Uint8Array> {
  const parts: Uint8Array[] = [];
  for (let i = 0; i < bytes.length; i += blockInput) {
    parts.push(await bgzfBlock(bytes.slice(i, i + blockInput)));
  }
  parts.push(BGZF_EOF);
  return new Uint8Array(Buffer.concat(parts));
}

describe('sniffCompression', () => {
  it('tells plain text, gzip and BGZF apart', async () => {
    expect(sniffCompression(SOURCE)).toBe('none');
    expect(sniffCompression(new Uint8Array(gzipSync(SOURCE)))).toBe('gzip');
    expect(sniffCompression(await bgzipped(SOURCE))).toBe('bgzf');
  });
});

describe('bgzfMembers', () => {
  it('splits a BGZF stream into whole members, whatever the chunking', async () => {
    const file = await bgzipped(SOURCE, 300);
    const sizes: number[] = [];
    for await (const member of bgzfMembers(chunked(file, 7))) {
      expect(bgzfBlockSize(member, 0)).toBe(member.length);
      sizes.push(member.length);
    }
    expect(sizes.length).toBeGreaterThan(3); // several data blocks plus the EOF block
    expect(sizes.reduce((a, b) => a + b, 0)).toBe(file.length);
  });

  it('refuses a truncated file', async () => {
    const file = await bgzipped(SOURCE, 300);
    const read = async () => {
      for await (const member of bgzfMembers(
        chunked(file.slice(0, file.length - 40), 64),
      )) {
        void member;
      }
    };
    await expect(read()).rejects.toThrow(/ends partway/);
  });
});

describe('digestFastaChunks', () => {
  it('matches the refget CLI digest for a plain FASTA, in any chunk size', async () => {
    for (const size of [1, 13, 4096]) {
      const result = await digestFastaChunks(chunked(SOURCE, size), wasm);
      expect(result.digest).toBe(DIGEST);
    }
  });

  it('reports the sequences in file order', async () => {
    const result = await digestFastaChunks(chunked(SOURCE, 512), wasm);
    expect(result.n_sequences).toBe(5);
    expect(result.sequences.map((s) => s.name)).toEqual([
      'chrA',
      'chrB',
      'chrC',
      'chrD',
      'chrE',
    ]);
    expect(result.sequences[0].sha512t24u).toBe('-YdTIx-7BHGG-SwIgibnGX3Ze5ckvvsW');
  });

  it('reads gzip and bgzip files to the same digest', async () => {
    const gz = new Uint8Array(gzipSync(SOURCE));
    expect((await digestFastaChunks(chunked(gz, 50), wasm)).digest).toBe(DIGEST);
    const bgz = await bgzipped(SOURCE, 300);
    expect((await digestFastaChunks(chunked(bgz, 50), wasm)).digest).toBe(DIGEST);
  });

  it('reports bytes read', async () => {
    const seen: number[] = [];
    await digestFastaChunks(chunked(SOURCE, 500), wasm, {
      onProgress: (n) => seen.push(n),
    });
    expect(seen[seen.length - 1]).toBe(SOURCE.length);
  });

  it('rejects empty and non-FASTA input', async () => {
    await expect(
      digestFastaChunks(chunked(new Uint8Array(0), 10), wasm),
    ).rejects.toThrow(/empty/);
    const text = new TextEncoder().encode('hello world\nnot a fasta\n');
    await expect(digestFastaChunks(chunked(text, 10), wasm)).rejects.toThrow();
  });
});

describe('looksLikeFasta', () => {
  it('accepts FASTA names, plain or gzipped', () => {
    for (const name of [
      'a.fa',
      'b.fasta',
      'c.fna',
      'd.fa.gz',
      'E.FASTA.GZ',
      'f.fa.bgz',
    ]) {
      expect(looksLikeFasta(name)).toBe(true);
    }
    expect(looksLikeFasta('notes.txt')).toBe(false);
    expect(looksLikeFasta('reads.bam')).toBe(false);
  });
});
