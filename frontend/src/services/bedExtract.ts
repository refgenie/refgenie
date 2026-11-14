/**
 * Read many BED regions of a genome in the browser, a few at a time, with the
 * same byte-range reader as the single-region box.
 */
import { extractRegion, type RegionSeq } from './regionExtract';
import { bedFastaRecord, type BedRegion } from '../utils/parseBed';

/** Reads in flight at once: enough to hide latency, few enough to stay polite. */
export const BED_CONCURRENCY = 6;

export interface ExtractBedArgs {
  baseUrl: string;
  regions: readonly BedRegion[];
  /** The genome's sequences by their own name. */
  sequences: ReadonlyMap<string, RegionSeq>;
  concurrency?: number;
  signal?: AbortSignal;
  /** Called after each region, with the count done so far. */
  onProgress?: (done: number, total: number) => void;
  /** For tests: the single-region reader. */
  read?: typeof extractRegion;
}

/** One FASTA record per region, in BED order. */
export async function extractBed({
  baseUrl,
  regions,
  sequences,
  concurrency = BED_CONCURRENCY,
  signal,
  onProgress,
  read = extractRegion,
}: ExtractBedArgs): Promise<string[]> {
  const records = new Array<string>(regions.length);
  // One failure stops the other readers too.
  const stop = new AbortController();
  const forward = () => stop.abort(signal?.reason);
  if (signal?.aborted) forward();
  signal?.addEventListener('abort', forward);
  const inner = stop.signal;
  let next = 0;
  let done = 0;
  const worker = async () => {
    while (next < regions.length) {
      inner.throwIfAborted();
      const index = next;
      next += 1;
      const region = regions[index];
      const seq = sequences.get(region.name);
      if (!seq)
        throw new Error(`Line ${region.line}: no sequence named "${region.chrom}".`);
      let bases: string;
      try {
        bases = await read({
          baseUrl,
          seq,
          start: region.start,
          end: region.end,
          signal: inner,
        });
      } catch (error) {
        if (inner.aborted) throw error;
        stop.abort();
        const message = error instanceof Error ? error.message : String(error);
        throw new Error(`Line ${region.line} (${region.header}): ${message}`);
      }
      records[index] = bedFastaRecord(region.header, bases);
      done += 1;
      onProgress?.(done, regions.length);
    }
  };
  try {
    const workers = Array.from(
      { length: Math.min(concurrency, regions.length) },
      worker,
    );
    await Promise.all(workers);
  } catch (error) {
    stop.abort();
    throw error;
  } finally {
    signal?.removeEventListener('abort', forward);
  }
  signal?.throwIfAborted();
  return records;
}
