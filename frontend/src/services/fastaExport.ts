/**
 * Rebuild a genome's FASTA in the browser from its refget store.
 *
 * An Encoded refget store keeps each sequence as a headerless bit-packed
 * `.seq` file, so a window of bases maps to a byte range. This module reads
 * those byte ranges over HTTP `Range:` requests and decodes them with the
 * published `@databio/gtars` wasm primitives, one window at a time, so memory
 * stays at about one window whatever the genome size.
 *
 * The output matches `refget store export -l name`: `>name` headers (no
 * description), uppercase bases, 80-column lines.
 *
 * No React or DOM here, so it runs under vitest in Node. The browser sink is
 * `saveStream.ts`.
 */

import type { SequenceResponse } from '../types/api';

type Gtars = typeof import('@databio/gtars');
export type FastaWasm = Pick<Gtars, 'encodedByteRange' | 'decodeEncodedRange'>;

/** One row of a collection `.rgsi`, in collection (= FASTA) order. */
export interface FastaSeq {
  name: string;
  length: number;
  alphabet: string;
  description: string | null;
  sha512t24u: string;
}

export interface FastaCollection {
  seqTemplate: string;
  seqs: FastaSeq[];
}

export const DEFAULT_SEQ_TEMPLATE = 'sequences/%s2/%s.seq';
export const DEFAULT_COLLECTION_TEMPLATE = 'collections/%s.rgsi';
const DEFAULT_LINE_WIDTH = 80;
const DEFAULT_WINDOW_BASES = 8_000_000;
const NEWLINE = 0x0a;

/** Alphabets whose `.seq` bytes are bit-packed and can be decoded by range. */
export const DECODABLE = new Set(['dna2bit', 'dna3bit', 'dnaio', 'dnaiupac']);

/**
 * Parse a collection `.rgsi`: `name, length, alphabet, sha512t24u, md5,
 * description`. Lines starting with `#` are headers.
 */
export function parseRgsiFull(text: string): FastaSeq[] {
  const seqs: FastaSeq[] = [];
  for (const line of text.split('\n')) {
    const row = line.replace(/\r$/, '');
    if (row === '' || row.startsWith('#')) continue;
    const cols = row.split('\t');
    const description = cols.slice(5).join('\t');
    seqs.push({
      name: cols[0],
      length: Number(cols[1]),
      alphabet: cols[2] ?? '',
      sha512t24u: cols[3] ?? '',
      description: description === '' ? null : description,
    });
  }
  return seqs;
}

/** Expand a refget store path template, the same rules as gtars `expand_template`. */
export function expandTemplate(template: string, digest: string): string {
  return template
    .split('%s2')
    .join(digest.slice(0, 2))
    .split('%s4')
    .join(digest.slice(0, 4))
    .split('%s')
    .join(digest);
}

function withSlash(baseUrl: string): string {
  return baseUrl.endsWith('/') ? baseUrl : `${baseUrl}/`;
}

async function fetchOk(url: string, fetchImpl: typeof fetch, init?: RequestInit) {
  const response = await fetchImpl(url, init);
  if (!response.ok) throw new Error(`Could not fetch ${url} (HTTP ${response.status})`);
  return response;
}

/**
 * The FASTA sequence list from the API's `/sequences` rows, or null if any row
 * lacks an alphabet (then the store's `.rgsi` is read instead).
 */
export function fastaSeqsFromRows(
  rows: Pick<
    SequenceResponse,
    'name' | 'length' | 'alphabet' | 'description' | 'sha512t24u'
  >[],
): FastaSeq[] | null {
  const seqs: FastaSeq[] = [];
  for (const row of rows) {
    if (!row.alphabet) return null;
    seqs.push({
      name: row.name,
      length: row.length,
      alphabet: row.alphabet,
      description: row.description || null,
      sha512t24u: row.sha512t24u,
    });
  }
  return seqs;
}

/**
 * Read the store manifest and the collection's sequence list. Pass `known`
 * (e.g. from {@link fastaSeqsFromRows}) to skip fetching the `.rgsi`.
 */
export async function loadCollection(
  baseUrl: string,
  digest: string,
  fetchImpl: typeof fetch = fetch,
  known?: FastaSeq[] | null,
): Promise<FastaCollection> {
  const base = withSlash(baseUrl);
  const manifest = (await (
    await fetchOk(new URL('rgstore.json', base).href, fetchImpl)
  ).json()) as { seqdata_path_template?: string; collections_path_template?: string };
  const seqTemplate = manifest.seqdata_path_template || DEFAULT_SEQ_TEMPLATE;
  if (known) return { seqTemplate, seqs: known };
  const collTemplate =
    manifest.collections_path_template || DEFAULT_COLLECTION_TEMPLATE;
  const collUrl = new URL(expandTemplate(collTemplate, digest), base).href;
  const text = await (await fetchOk(collUrl, fetchImpl)).text();
  return { seqTemplate, seqs: parseRgsiFull(text) };
}

/** Throw if any sequence uses an alphabet that cannot be decoded by range. */
/** True if the browser can rebuild this sequence's bases. */
export function isDecodable(seq: Pick<FastaSeq, 'alphabet'>): boolean {
  return DECODABLE.has(seq.alphabet.toLowerCase());
}

export function assertDecodable(seqs: FastaSeq[]): void {
  const bad = seqs.find((seq) => !isDecodable(seq));
  if (bad) {
    throw new Error(
      `Sequence "${bad.name}" uses the ${bad.alphabet || 'unknown'} alphabet, ` +
        'which cannot be rebuilt as FASTA in the browser.',
    );
  }
}

const encoder = new TextEncoder();

function headerLine(seq: FastaSeq): string {
  return `>${seq.name}\n`;
}

/** The exact byte count `fastaChunks` will produce. */
export function fastaSize(seqs: FastaSeq[], lineWidth = DEFAULT_LINE_WIDTH): number {
  let total = 0;
  for (const seq of seqs) {
    total += encoder.encode(headerLine(seq)).length;
    total += seq.length + Math.ceil(seq.length / lineWidth);
  }
  return total;
}

interface Window {
  seq: FastaSeq;
  url: string;
  start: number;
  end: number;
  /** True for the first window of a sequence. */
  first: boolean;
  /** True for the last window of a sequence. */
  last: boolean;
}

export interface FastaChunksOptions {
  baseUrl: string;
  digest: string;
  wasm: FastaWasm;
  fetchImpl?: typeof fetch;
  lineWidth?: number;
  windowBases?: number;
  onProgress?: (basesDone: number, basesTotal: number) => void;
  signal?: AbortSignal;
  /** An already loaded collection, to skip refetching it. */
  collection?: FastaCollection;
}

/** The `.seq` URL of one sequence in a store. */
export function seqUrl(baseUrl: string, seqTemplate: string, digest: string): string {
  return new URL(expandTemplate(seqTemplate, digest), withSlash(baseUrl)).href;
}

/**
 * Read bases [start, end) (0-based, half-open) of one encoded `.seq` file:
 * fetch only the covering bytes with an HTTP `Range:` request and decode them.
 * Fails if the server ignores `Range`, rather than fetching whole files per window.
 */
export async function readEncodedBases({
  url,
  start,
  end,
  alphabet,
  wasm,
  fetchImpl = fetch,
  signal,
}: {
  url: string;
  start: number;
  end: number;
  alphabet: string;
  wasm: FastaWasm;
  fetchImpl?: typeof fetch;
  signal?: AbortSignal;
}): Promise<string> {
  if (end <= start) return '';
  const [bs, be] = wasm.encodedByteRange(start, end, alphabet);
  const response = await fetchImpl(url, {
    headers: { Range: `bytes=${bs}-${be - 1}` },
    signal,
  });
  let bytes: Uint8Array;
  if (response.status === 206) {
    bytes = new Uint8Array(await response.arrayBuffer());
    if (bytes.length < be - bs) throw new Error(`Short read from ${url}`);
  } else if (response.status === 200) {
    // Reading on would download the whole file once per window.
    await response.body?.cancel().catch(() => {});
    throw new Error(
      'The store server ignored a byte-range request, so sequences cannot be ' +
        `read from it piece by piece (${url}).`,
    );
  } else {
    throw new Error(`Could not fetch ${url} (HTTP ${response.status})`);
  }
  return wasm.decodeEncodedRange(bytes, bs, start, end, alphabet);
}

/**
 * Stream a genome's FASTA as UTF-8 chunks: one header chunk per sequence, then
 * one chunk per decoded window. One window fetch stays in flight ahead of the
 * window being decoded, so network and decode overlap.
 */
export async function* fastaChunks(
  opts: FastaChunksOptions,
): AsyncGenerator<Uint8Array> {
  const {
    baseUrl,
    digest,
    wasm,
    fetchImpl = fetch,
    lineWidth = DEFAULT_LINE_WIDTH,
    windowBases = DEFAULT_WINDOW_BASES,
    onProgress,
    signal,
  } = opts;
  if (lineWidth < 1) throw new Error('lineWidth must be at least 1');
  if (windowBases < 1) throw new Error('windowBases must be at least 1');

  const base = withSlash(baseUrl);
  const { seqTemplate, seqs } =
    opts.collection ?? (await loadCollection(base, digest, fetchImpl));
  assertDecodable(seqs);
  signal?.throwIfAborted();

  const basesTotal = seqs.reduce((sum, seq) => sum + seq.length, 0);
  let basesDone = 0;

  // Sequences with no bases still need their header, so they get one empty window.
  const windows: Window[] = [];
  for (const seq of seqs) {
    const url = seqUrl(base, seqTemplate, seq.sha512t24u);
    if (seq.length === 0) {
      windows.push({ seq, url, start: 0, end: 0, first: true, last: true });
      continue;
    }
    for (let start = 0; start < seq.length; start += windowBases) {
      const end = Math.min(start + windowBases, seq.length);
      windows.push({
        seq,
        url,
        start,
        end,
        first: start === 0,
        last: end === seq.length,
      });
    }
  }

  const fetchWindow = (win: Window): Promise<string> =>
    win.end === win.start
      ? Promise.resolve('')
      : readEncodedBases({
          url: win.url,
          start: win.start,
          end: win.end,
          alphabet: win.seq.alphabet,
          wasm,
          fetchImpl,
          signal,
        });

  const prefetch = (index: number): Promise<string> | null => {
    if (index >= windows.length) return null;
    const pending = fetchWindow(windows[index]);
    // Awaited later; this only stops an abort from surfacing as unhandled.
    pending.catch(() => {});
    return pending;
  };

  let column = 0;
  let next = prefetch(0);
  for (let i = 0; i < windows.length; i += 1) {
    signal?.throwIfAborted();
    const win = windows[i];
    const current = next as Promise<string>;
    next = prefetch(i + 1);

    if (win.first) {
      column = 0;
      yield encoder.encode(headerLine(win.seq));
    }

    const bases = await current;
    signal?.throwIfAborted();
    const out = wrapBases(bases, column, lineWidth, win.last);
    column = out.column;
    if (out.bytes.length > 0) yield out.bytes;

    basesDone += win.end - win.start;
    onProgress?.(basesDone, basesTotal);
  }
}

/**
 * Encode a window of bases with a newline every `lineWidth` bases, carrying the
 * column across windows. A sequence's last window ends its final partial line.
 */
function wrapBases(
  bases: string,
  column: number,
  lineWidth: number,
  last: boolean,
): { bytes: Uint8Array; column: number } {
  const n = bases.length;
  const breaks = Math.floor((column + n) / lineWidth);
  const tail = last && (column + n) % lineWidth !== 0 ? 1 : 0;
  const bytes = new Uint8Array(n + breaks + tail);
  let pos = 0;
  let i = 0;
  while (i < n) {
    const take = Math.min(lineWidth - column, n - i);
    encoder.encodeInto(bases.slice(i, i + take), bytes.subarray(pos, pos + take));
    pos += take;
    i += take;
    column += take;
    if (column === lineWidth) {
      bytes[pos++] = NEWLINE;
      column = 0;
    }
  }
  if (tail) {
    bytes[pos++] = NEWLINE;
    column = 0;
  }
  return { bytes, column };
}

let gtars: Gtars | null = null;

/** Load and initialise the gtars wasm module once, on first use. */
export async function loadGtars(): Promise<Gtars> {
  if (gtars) return gtars;
  const mod = await import('@databio/gtars');
  await mod.default();
  gtars = mod;
  return gtars;
}
