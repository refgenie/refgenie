import { gzipSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import {
  bedFastaFileName,
  bedFastaRecord,
  bedHeader,
  bedOverCap,
  exampleBed,
  parseBed,
  readBedFile,
} from './parseBed';

const SEQS = [
  {
    name: 'chr1',
    length: 2_000_000,
    aliases: [{ alias: 'NC_000001.11' }],
    alphabet: 'dna2bit',
  },
  { name: 'chrM', length: 16_569, aliases: [{ alias: 'MT' }], alphabet: 'dna2bit' },
  { name: 'chrUn', length: 500, alphabet: 'protein' },
];
const readable = (s: { alphabet?: string | null }) => s.alphabet === 'dna2bit';

describe('bedHeader', () => {
  it('writes chr:start-end, or name::chr:start-end, in BED coordinates', () => {
    expect(bedHeader('chrM', 0, 60)).toBe('chrM:0-60');
    expect(bedHeader('chr1', 999_999, 1_000_500, 'peak1')).toBe(
      'peak1::chr1:999999-1000500',
    );
  });

  it('makes one-line FASTA records', () => {
    expect(bedFastaRecord('chrM:0-4', 'GATC')).toBe('>chrM:0-4\nGATC\n');
  });
});

describe('parseBed', () => {
  it('reads tabs or spaces, skips headers, and keeps the name column', () => {
    const bed = parseBed(
      [
        'track name=peaks',
        'browser position chr1:1-100',
        '# a comment',
        '',
        'chrM\t0\t60',
        'chr1 999999 1000500 peak1 0 +',
        'chr1\t10\t20\t.\r',
      ].join('\n'),
      SEQS,
      readable,
    );
    expect(bed.problems).toEqual([]);
    expect(bed.skipped).toBe(4);
    expect(bed.regions.map((r) => [r.line, r.header])).toEqual([
      [5, 'chrM:0-60'],
      [6, 'peak1::chr1:999999-1000500'],
      [7, 'chr1:10-20'],
    ]);
    expect(bed.bases).toBe(60 + 501 + 10);
  });

  it('matches aliases, keeping the name as written for the header', () => {
    const bed = parseBed('MT\t0\t10\nNC_000001.11\t5\t9\n', SEQS, readable);
    expect(bed.regions.map((r) => [r.name, r.header])).toEqual([
      ['chrM', 'MT:0-10'],
      ['chr1', 'NC_000001.11:5-9'],
    ]);
  });

  it('keeps a tab-separated name with spaces', () => {
    const bed = parseBed('chrM\t0\t10\tmy peak\n', SEQS, readable);
    expect(bed.regions[0].header).toBe('my peak::chrM:0-10');
  });

  it('lists problem lines with their line numbers and lets the rest run', () => {
    const bed = parseBed(
      [
        'chrM\t0\t60',
        'chrm\t0\t10',
        'chr1\t50\t50',
        'chrM\t16000\t17000',
        'chr1\tabc\t10',
        'chr1\t10',
        'chrUn\t0\t10',
        'chr2\t0\t10',
      ].join('\n'),
      SEQS,
      readable,
    );
    expect(bed.regions).toHaveLength(1);
    expect(bed.problems.map((p) => [p.line, p.message])).toEqual([
      [2, 'No sequence named "chrm". Did you mean "chrM"?'],
      [3, 'The start (50) must be less than the end (50).'],
      [4, 'The end (17,000) is past the end of chrM (16,569 bases).'],
      [5, 'Not a BED line: it needs a sequence name, a start and an end.'],
      [6, 'Not a BED line: it needs a sequence name, a start and an end.'],
      [7, 'chrUn uses the protein alphabet, which cannot be read here.'],
      [8, 'No sequence named "chr2".'],
    ]);
    expect(bed.problems[0].text).toBe('chrm\t0\t10');
  });
});

describe('bedOverCap', () => {
  it('says which limit a file passes', () => {
    expect(bedOverCap({ regions: new Array(3), bases: 10 }, 5, 100)).toBeNull();
    expect(bedOverCap({ regions: new Array(6), bases: 10 }, 5, 100)).toBe(
      'This file has 6 regions (the limit is 5).',
    );
    expect(bedOverCap({ regions: new Array(6), bases: 200 }, 5, 100)).toBe(
      'This file has 6 regions (the limit is 5) and 200 bases (the limit is 100).',
    );
  });
});

describe('exampleBed', () => {
  it('builds valid regions from the genome itself', () => {
    const text = exampleBed(SEQS) as string;
    expect(text.split('\n').filter(Boolean)).toEqual([
      'chr1\t999999\t1000500\tregion_1',
      'chrM\t8284\t8584\tregion_2',
      'chrUn\t250\t500\tregion_3',
      'chrUn\t0\t60',
    ]);
    const bed = parseBed(text, SEQS);
    expect(bed.problems).toEqual([]);
    expect(bed.regions).toHaveLength(4);
  });

  it('works for a genome of one sequence', () => {
    expect(exampleBed([{ name: 'x', length: 100 }])).toBe('x\t0\t100\tregion_1\n');
    expect(exampleBed([])).toBeNull();
  });
});

describe('bedFastaFileName', () => {
  it('names the FASTA after the BED file', () => {
    expect(bedFastaFileName('peaks.bed.gz')).toBe('peaks.fa');
    expect(bedFastaFileName('my peaks.bed')).toBe('my_peaks.fa');
    expect(bedFastaFileName(null)).toBe('regions.fa');
  });
});

describe('readBedFile', () => {
  it('reads plain text', async () => {
    expect(await readBedFile(new Blob(['chrM\t0\t10\n']))).toBe('chrM\t0\t10\n');
  });

  it('unzips gzip by its magic bytes', async () => {
    const zipped = new Blob([gzipSync('chrM\t0\t10\n')]);
    expect(await readBedFile(zipped)).toBe('chrM\t0\t10\n');
  });
});
