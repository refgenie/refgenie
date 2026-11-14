/**
 * Read one region of a genome in the browser from its refget store.
 *
 * Only the bytes covering the region are fetched (an HTTP `Range:` request on
 * the sequence's bit-packed `.seq` file) and decoded with the same
 * `@databio/gtars` wasm module the FASTA download uses, so a 500-base read
 * costs a few hundred bytes, not a chromosome.
 *
 * Ported from the RefgetStore Explorer's extractSequence.ts, which wraps
 * gtars' `RemoteRefgetStore`. That class caches in OPFS through
 * `navigator.storage`, which browsers only offer on https or localhost; a
 * refgenie server on plain http would fail. This module reuses the FASTA
 * download's range reader instead: same wasm calls, one wasm copy, and a
 * `fetch` that tests can swap.
 */

import {
  DECODABLE,
  loadGtars,
  loadSeqTemplate,
  readEncodedBases,
  seqUrl,
  type FastaWasm,
} from './fastaExport';

export interface RegionSeq {
  name: string;
  length: number;
  alphabet: string | null;
  sha512t24u: string;
}

/** True if a sequence's bases can be read by range in the browser. */
export function isExtractable(alphabet: string | null | undefined): boolean {
  return !!alphabet && DECODABLE.has(alphabet.toLowerCase());
}

export interface ExtractRegionArgs {
  baseUrl: string;
  seq: RegionSeq;
  /** 0-based, inclusive. */
  start: number;
  /** 0-based, exclusive. */
  end: number;
  wasm?: FastaWasm;
  fetchImpl?: typeof fetch;
  signal?: AbortSignal;
}

/** The bases [start, end) of `seq`, uppercase as stored. */
export async function extractRegion({
  baseUrl,
  seq,
  start,
  end,
  wasm,
  fetchImpl = fetch,
  signal,
}: ExtractRegionArgs): Promise<string> {
  if (!isExtractable(seq.alphabet)) {
    throw new Error(
      `${seq.name} uses the ${seq.alphabet ?? 'unknown'} alphabet, which cannot be read here.`,
    );
  }
  if (start < 0 || end > seq.length || start > end) {
    throw new Error(
      `Region ${start}-${end} is outside ${seq.name} (${seq.length} bases).`,
    );
  }
  const [template, gtars] = await Promise.all([
    loadSeqTemplate(baseUrl, fetchImpl),
    wasm ?? loadGtars(),
  ]);
  return readEncodedBases({
    url: seqUrl(baseUrl, template, seq.sha512t24u),
    start,
    end,
    alphabet: seq.alphabet as string,
    wasm: gtars,
    fetchImpl,
    signal,
  });
}
