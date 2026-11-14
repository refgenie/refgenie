/**
 * Read one region of a genome in the browser from its refget store.
 *
 * Only the bytes covering the region are fetched (an HTTP `Range:` request on
 * the sequence's bit-packed `.seq` file) and decoded with the same
 * `@databio/gtars` wasm module the FASTA download uses, so a 500-base read
 * costs a few hundred bytes, not a chromosome. Not gtars' `RemoteRefgetStore`:
 * it caches in OPFS, which browsers offer only on https or localhost.
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
