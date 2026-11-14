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
  isDecodable,
  loadGtars,
  loadSeqTemplate,
  readEncodedBases,
  seqUrl,
  undecodableMessage,
  type FastaWasm,
} from './fastaExport';

export interface RegionSeq {
  name: string;
  length: number;
  alphabet: string | null;
  sha512t24u: string;
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
  if (!isDecodable(seq.alphabet)) throw new Error(undecodableMessage(seq));
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
