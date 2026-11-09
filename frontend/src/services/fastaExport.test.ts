// @vitest-environment node
/// <reference types="node" />
import { readFileSync, existsSync } from 'node:fs';
import { createRequire } from 'node:module';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { beforeAll, describe, expect, it } from 'vitest';
import {
  assertDecodable,
  fastaChunks,
  fastaSeqsFromRows,
  fastaSize,
  loadCollection,
  parseRgsiFull,
  type FastaWasm,
} from './fastaExport';
import { extractRegion } from './regionExtract';

const FIXTURE = join(
  dirname(fileURLToPath(import.meta.url)),
  '../test/fixtures/fasta-store',
);
const STORE = join(FIXTURE, 'store');
const BASE = 'https://store.test/s/';
// The browser writes `>name` headers, like `refget store export -l name`, so
// the fixture's header descriptions are dropped.
const SOURCE = Buffer.from(
  readFileSync(join(FIXTURE, 'source.fa'), 'utf8').replace(/^(>\S+)[^\n]*$/gm, '$1'),
);
const DIGEST = readFileSync(join(FIXTURE, 'digest.txt'), 'utf8').trim();

let wasm: FastaWasm;

beforeAll(async () => {
  // The package `exports` map does not expose the .wasm file, so resolve the
  // main entry and read the binary next to it.
  const dir = dirname(createRequire(import.meta.url).resolve('@databio/gtars'));
  const mod = await import('@databio/gtars');
  await mod.default({ module_or_path: readFileSync(join(dir, 'gtars_js_bg.wasm')) });
  wasm = mod;
});

function storeFile(url: string): Buffer | null {
  const path = new URL(url).pathname;
  const prefix = new URL(BASE).pathname;
  if (!path.startsWith(prefix)) return null;
  const file = join(STORE, decodeURIComponent(path.slice(prefix.length)));
  return existsSync(file) ? readFileSync(file) : null;
}

function requestUrl(input: RequestInfo | URL): string {
  return typeof input === 'string'
    ? input
    : input instanceof URL
      ? input.href
      : input.url;
}

/** Serves the fixture store under BASE, honouring `Range: bytes=a-b`. */
const fixtureFetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
  const body = storeFile(requestUrl(input));
  if (!body) return new Response('missing', { status: 404 });
  const range = new Headers(init?.headers).get('Range');
  const match = range?.match(/^bytes=(\d+)-(\d+)$/);
  if (match) {
    const slice = body.subarray(Number(match[1]), Number(match[2]) + 1);
    return new Response(new Uint8Array(slice), { status: 206 });
  }
  return new Response(new Uint8Array(body), { status: 200 });
}) as typeof fetch;

/** Like fixtureFetch, but ignores Range and always sends the whole file. */
const noRangeFetch = (async (input: RequestInfo | URL) => {
  const body = storeFile(requestUrl(input));
  if (!body) return new Response('missing', { status: 404 });
  return new Response(new Uint8Array(body), { status: 200 });
}) as typeof fetch;

async function collect(gen: AsyncGenerator<Uint8Array>): Promise<Buffer> {
  const parts: Uint8Array[] = [];
  for await (const chunk of gen) parts.push(chunk);
  return Buffer.concat(parts);
}

describe('fastaChunks', () => {
  it('writes name-only headers', async () => {
    const out = await collect(
      fastaChunks({ baseUrl: BASE, digest: DIGEST, wasm, fetchImpl: fixtureFetch }),
    );
    expect(out.toString('utf8')).toContain('>chrB\n');
    expect(out.toString('utf8')).not.toContain('test contig B');
  });

  it('exports only the sequences in the collection it is given', async () => {
    const full = await loadCollection(BASE, DIGEST, fixtureFetch);
    const collection = {
      ...full,
      seqs: full.seqs.filter((s) => ['chrB', 'chrE'].includes(s.name)),
    };
    const out = await collect(
      fastaChunks({
        baseUrl: BASE,
        digest: DIGEST,
        wasm,
        fetchImpl: fixtureFetch,
        collection,
      }),
    );
    const headers = out.toString('utf8').match(/^>.*$/gm);
    expect(headers).toEqual(['>chrB', '>chrE']);
    expect(out.length).toBe(fastaSize(collection.seqs));
  });

  it('rebuilds the source FASTA byte for byte', async () => {
    const out = await collect(
      fastaChunks({ baseUrl: BASE, digest: DIGEST, wasm, fetchImpl: fixtureFetch }),
    );
    expect(out.equals(SOURCE)).toBe(true);
  });

  it.each([1, 3, 37, 79, 80, 81])(
    'is byte-identical with %i-base windows',
    async (windowBases) => {
      const out = await collect(
        fastaChunks({
          baseUrl: BASE,
          digest: DIGEST,
          wasm,
          fetchImpl: fixtureFetch,
          windowBases,
        }),
      );
      expect(out.toString('utf8')).toBe(SOURCE.toString('utf8'));
    },
  );

  it('fails clearly when the server ignores Range', async () => {
    await expect(
      collect(
        fastaChunks({
          baseUrl: BASE,
          digest: DIGEST,
          wasm,
          fetchImpl: noRangeFetch,
          windowBases: 37,
        }),
      ),
    ).rejects.toThrow(/ignored a byte-range request/);
  });

  it('reports progress up to the total base count', async () => {
    const seen: Array<[number, number]> = [];
    await collect(
      fastaChunks({
        baseUrl: BASE,
        digest: DIGEST,
        wasm,
        fetchImpl: fixtureFetch,
        windowBases: 100,
        onProgress: (done, total) => seen.push([done, total]),
      }),
    );
    expect(seen[seen.length - 1]).toEqual([1451, 1451]);
  });

  it('stops with an AbortError when the signal is aborted', async () => {
    const controller = new AbortController();
    const gen = fastaChunks({
      baseUrl: BASE,
      digest: DIGEST,
      wasm,
      fetchImpl: fixtureFetch,
      windowBases: 10,
      signal: controller.signal,
    });
    await gen.next();
    controller.abort();
    await expect(collect(gen)).rejects.toMatchObject({ name: 'AbortError' });
  });
});

describe('fastaSize', () => {
  it('equals the source FASTA length', async () => {
    const { seqs } = await loadCollection(BASE, DIGEST, fixtureFetch);
    expect(fastaSize(seqs)).toBe(SOURCE.length);
  });
});

describe('fastaSeqsFromRows', () => {
  it('lets loadCollection skip the .rgsi', async () => {
    const fromRgsi = await loadCollection(BASE, DIGEST, fixtureFetch);
    const known = fastaSeqsFromRows(fromRgsi.seqs);
    const urls: string[] = [];
    const spy: typeof fetch = (input, init) => {
      urls.push(String(input));
      return fixtureFetch(input, init);
    };
    const { seqs } = await loadCollection(BASE, DIGEST, spy, known);
    expect(seqs).toEqual(fromRgsi.seqs);
    expect(urls.some((u) => u.endsWith('.rgsi'))).toBe(false);
  });

  it('gives null when a row has no alphabet', () => {
    const row = {
      name: 'c',
      length: 1,
      alphabet: null,
      description: null,
      sha512t24u: 'x',
    };
    expect(fastaSeqsFromRows([row])).toBeNull();
  });
});

describe('parseRgsiFull and assertDecodable', () => {
  it('keeps file order and reads descriptions', async () => {
    const { seqs } = await loadCollection(BASE, DIGEST, fixtureFetch);
    expect(seqs.map((s) => s.name)).toEqual(['chrA', 'chrB', 'chrC', 'chrD', 'chrE']);
    expect(seqs[0].description).toBeNull();
    expect(seqs[1].description).toBe('test contig B');
  });

  it('names the first sequence with an unsupported alphabet', () => {
    const seqs = parseRgsiFull(
      '#name\tlength\talphabet\tsha512t24u\tmd5\tdescription\n' +
        'chr1\t10\tdna2bit\tAAAA\tmd5\t\n' +
        'prot1\t5\tprotein\tBBBB\tmd5\t\n',
    );
    expect(() => assertDecodable(seqs)).toThrow(/prot1.*protein/);
  });
});

describe('extractRegion', () => {
  // Bases of each source record, for comparison.
  const sourceBases = new Map<string, string>();
  for (const record of SOURCE.toString('utf8').split('>').slice(1)) {
    const [name, ...lines] = record.split('\n');
    sourceBases.set(name, lines.join(''));
  }

  it('reads exact bases for regions of every sequence', async () => {
    const { seqs } = await loadCollection(BASE, DIGEST, fixtureFetch);
    for (const seq of seqs) {
      const bases = sourceBases.get(seq.name) as string;
      const spans: [number, number][] = [
        [0, seq.length],
        [0, 1],
        [seq.length - 1, seq.length],
        [Math.floor(seq.length / 3), Math.floor((2 * seq.length) / 3)],
        [1, Math.min(seq.length, 82)],
      ];
      for (const [start, end] of spans) {
        const out = await extractRegion({
          baseUrl: BASE,
          seq,
          start,
          end,
          wasm,
          fetchImpl: fixtureFetch,
        });
        expect(out, `${seq.name}:${start}-${end}`).toBe(bases.slice(start, end));
      }
    }
  });

  it('fetches only the bytes the region needs', async () => {
    const { seqs } = await loadCollection(BASE, DIGEST, fixtureFetch);
    const seq = seqs.find((s) => s.name === 'chrD') as (typeof seqs)[number];
    const ranges: string[] = [];
    const spy = (async (input: RequestInfo | URL, init?: RequestInit) => {
      const range = new Headers(init?.headers).get('Range');
      if (range) ranges.push(range);
      return fixtureFetch(input, init);
    }) as typeof fetch;
    const out = await extractRegion({
      baseUrl: BASE,
      seq,
      start: 400,
      end: 420,
      wasm,
      fetchImpl: spy,
    });
    expect(out).toBe((sourceBases.get('chrD') as string).slice(400, 420));
    // dna2bit: 20 bases fit in 5-6 bytes.
    expect(ranges).toEqual(['bytes=100-104']);
  });

  it('refuses a range outside the sequence', async () => {
    const { seqs } = await loadCollection(BASE, DIGEST, fixtureFetch);
    await expect(
      extractRegion({
        baseUrl: BASE,
        seq: seqs[0],
        start: 0,
        end: seqs[0].length + 1,
        wasm,
        fetchImpl: fixtureFetch,
      }),
    ).rejects.toThrow(/outside/);
  });
});
