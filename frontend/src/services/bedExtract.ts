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
  /** Called after each batch, with the count done so far. */
  onProgress?: (done: number, total: number) => void;
  /** For tests: the single-region reader. */
  read?: typeof extractRegion;
}

/** One FASTA record per region, in BED order, reading `concurrency` at a time. */
export async function extractBed({
  baseUrl,
  regions,
  sequences,
  concurrency = BED_CONCURRENCY,
  signal,
  onProgress,
  read = extractRegion,
}: ExtractBedArgs): Promise<string[]> {
  const readOne = async (region: BedRegion) => {
    const seq = sequences.get(region.name);
    if (!seq) throw new Error(`Line ${region.line}: no sequence named "${region.chrom}".`);
    try {
      const { start, end } = region;
      return bedFastaRecord(region.header, await read({ baseUrl, seq, start, end, signal }));
    } catch (error) {
      if (signal?.aborted) throw error;
      const message = error instanceof Error ? error.message : String(error);
      throw new Error(`Line ${region.line} (${region.header}): ${message}`);
    }
  };
  const records: string[] = [];
  for (let i = 0; i < regions.length; i += concurrency) {
    signal?.throwIfAborted();
    records.push(...(await Promise.all(regions.slice(i, i + concurrency).map(readOne))));
    onProgress?.(records.length, regions.length);
  }
  return records;
}
