import { describe, expect, it, vi } from 'vitest';
import { extractBed } from './bedExtract';
import { parseBed } from '../utils/parseBed';
import type { RegionSeq } from './regionExtract';

const SEQS: RegionSeq[] = [
  { name: 'chrM', length: 100, alphabet: 'dna2bit', sha512t24u: 'm' },
  { name: 'chr1', length: 1000, alphabet: 'dna2bit', sha512t24u: 'one' },
];
const byName = new Map(SEQS.map((s) => [s.name, s]));
const fake = (start: number, end: number) => 'ACGT'.repeat(300).slice(start, end);

describe('extractBed', () => {
  it('reads every region, at most `concurrency` at once, in BED order', async () => {
    const bed = parseBed(
      ['chr1\t10\t20\tb', 'chrM\t0\t4', 'chr1\t0\t8', 'chrM\t4\t12', 'chr1\t5\t9'].join(
        '\n',
      ),
      SEQS,
    );
    let inFlight = 0;
    let most = 0;
    const read = vi.fn(async ({ start, end }: { start: number; end: number }) => {
      inFlight += 1;
      most = Math.max(most, inFlight);
      await new Promise((resolve) => setTimeout(resolve, 5 * (end % 3)));
      inFlight -= 1;
      return fake(start, end);
    });
    const progress: number[] = [];
    const records = await extractBed({
      baseUrl: 'https://store.test/',
      regions: bed.regions,
      sequences: byName,
      concurrency: 2,
      onProgress: (done) => progress.push(done),
      read: read as never,
    });
    expect(most).toBe(2);
    expect(progress).toEqual([1, 2, 3, 4, 5]);
    expect(records.join('')).toBe(
      '>b::chr1:10-20\nGTACGTACGT\n>chrM:0-4\nACGT\n>chr1:0-8\nACGTACGT\n' +
        '>chrM:4-12\nACGTACGT\n>chr1:5-9\nCGTA\n',
    );
  });

  it('names the line that failed, and stops the rest', async () => {
    const bed = parseBed('chrM\t0\t4\nchr1\t0\t4\nchr1\t4\t8\nchr1\t8\t12', SEQS);
    const read = vi.fn(async ({ seq }: { seq: RegionSeq }) => {
      if (seq.name === 'chr1') throw new Error('HTTP 500');
      return 'ACGT';
    });
    await expect(
      extractBed({
        baseUrl: 'https://store.test/',
        regions: bed.regions,
        sequences: byName,
        concurrency: 1,
        read: read as never,
      }),
    ).rejects.toThrow('Line 2 (chr1:0-4): HTTP 500');
    expect(read).toHaveBeenCalledTimes(2);
  });

  it('stops when cancelled', async () => {
    const bed = parseBed('chrM\t0\t4\nchrM\t4\t8\nchrM\t8\t12', SEQS);
    const abort = new AbortController();
    const read = vi.fn(async () => {
      abort.abort();
      return 'ACGT';
    });
    await expect(
      extractBed({
        baseUrl: 'https://store.test/',
        regions: bed.regions,
        sequences: byName,
        concurrency: 1,
        signal: abort.signal,
        read: read as never,
      }),
    ).rejects.toThrow();
    expect(read).toHaveBeenCalledTimes(1);
  });
});
