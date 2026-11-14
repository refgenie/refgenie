import { describe, expect, it } from 'vitest';
import {
  exampleRegion,
  MAX_EXTRACT_BASES,
  parseRegion,
  regionFasta,
  regionFileName,
  suggestSequences,
} from './parseRegion';

const SEQS = [
  {
    name: 'chr1',
    length: 248_956_422,
    aliases: [{ alias: 'NC_000001.11' }, { alias: '1' }],
  },
  { name: 'chr2', length: 242_193_529, aliases: [] },
  { name: 'chrM', length: 16_569 },
  { name: 'HLA-A*01:01:01:01', length: 3_503 },
];

describe('parseRegion', () => {
  it.each([
    ['chr1:1,000,000-1,000,500', 999_999, 1_000_500],
    ['chr1:1000000-1000500', 999_999, 1_000_500],
    [' chr1 : 1_000_000 - 1_000_500 ', 999_999, 1_000_500],
    ['chr1:1000000', 999_999, 1_000_000],
    ['chrM:1-60', 0, 60],
    ['chrM:16000-', 15_999, 16_569],
    ['chrM', 0, 16_569],
  ])('reads %s as 1-based, inclusive', (text, start, end) => {
    const region = parseRegion(text, SEQS);
    expect([region.start, region.end]).toEqual([start, end]);
  });

  it('labels the region in 1-based form', () => {
    expect(parseRegion('chr1:1,000,000-1,000,500', SEQS).label).toBe(
      'chr1:1000000-1000500',
    );
    expect(parseRegion('chrM', SEQS).label).toBe('chrM:1-16569');
  });

  it('resolves an alias to the sequence name', () => {
    const region = parseRegion('NC_000001.11:100-200', SEQS);
    expect(region.name).toBe('chr1');
    expect(region.label).toBe('chr1:100-200');
  });

  it('reads a name that contains colons', () => {
    expect(parseRegion('HLA-A*01:01:01:01', SEQS)).toMatchObject({
      start: 0,
      end: 3_503,
    });
    expect(parseRegion('HLA-A*01:01:01:01:10-20', SEQS)).toMatchObject({
      name: 'HLA-A*01:01:01:01',
      start: 9,
      end: 20,
    });
  });

  it('keeps spaces inside a name that takes coordinates', () => {
    const seqs = [{ name: 'contig one', length: 1_000 }];
    expect(parseRegion('contig one:10-20', seqs)).toMatchObject({
      name: 'contig one',
      start: 9,
      end: 20,
    });
    expect(parseRegion('contig one : 10 - 20', seqs)).toMatchObject({ start: 9, end: 20 });
    expect(() => parseRegion('contigone:10-20', seqs)).toThrow(/no sequence named/);
  });

  it('cuts an end past the sequence back to its length', () => {
    const region = parseRegion('chrM:16500-17000', SEQS);
    expect(region).toMatchObject({ start: 16_499, end: 16_569, clamped: true });
  });

  it('suggests close names for an unknown one', () => {
    expect(() => parseRegion('2:1-10', SEQS)).toThrow('Did you mean "chr2"?');
    expect(() => parseRegion('MT:1-10', SEQS)).toThrow('Did you mean "chrM"?');
    expect(() => parseRegion('CHR1:1-10', SEQS)).toThrow('Did you mean "chr1"?');
    expect(() => parseRegion('chrZ:1-10', SEQS)).toThrow(/no sequence named "chrZ"\.$/);
  });

  it('rejects a start past the end, a reversed range, and position 0', () => {
    expect(() => parseRegion('chrM:20000-20010', SEQS)).toThrow(/has 16,569 bases/);
    expect(() => parseRegion('chrM:200-100', SEQS)).toThrow(/before the start/);
    expect(() => parseRegion('chrM:0-10', SEQS)).toThrow(/start at 1/);
  });

  it('rejects text that is not a region', () => {
    expect(() => parseRegion('', SEQS)).toThrow(/Enter a region/);
    expect(() => parseRegion('chr1:abc', SEQS)).toThrow(/Cannot read/);
  });

  it('caps a region at the extraction limit', () => {
    expect(MAX_EXTRACT_BASES).toBe(100_000);
    expect(parseRegion(`chr1:1-${MAX_EXTRACT_BASES}`, SEQS).end).toBe(
      MAX_EXTRACT_BASES,
    );
    expect(() => parseRegion(`chr1:1-${MAX_EXTRACT_BASES + 1}`, SEQS)).toThrow(
      /in the browser is 100,000/,
    );
    expect(() => parseRegion('chr1', SEQS)).toThrow(/most you can extract/);
    // The code snippets lift the cap.
    expect(parseRegion('chr1', SEQS, Number.POSITIVE_INFINITY).start).toBe(0);
  });
});

describe('suggestSequences', () => {
  it('finds a name through an alias', () => {
    expect(suggestSequences('nc_000001.11', SEQS)).toEqual(['chr1']);
  });
});

describe('regionFasta', () => {
  it('writes a header and 80-column lines', () => {
    const bases = 'A'.repeat(80) + 'C'.repeat(5);
    expect(regionFasta('chr1:1-85', bases)).toBe(
      `>chr1:1-85\n${'A'.repeat(80)}\nCCCCC\n`,
    );
  });

  it('names the file after the region', () => {
    expect(regionFileName('chr1:1000000-1000500')).toBe('chr1_1000000-1000500.fa');
  });
});

describe('exampleRegion', () => {
  it('uses 1,000,000 on the longest sequence when it is long enough', () => {
    const seqs = [
      { name: 'small', length: 10 },
      { name: 'big', length: 5_000_000 },
    ];
    expect(exampleRegion(seqs)).toBe('big:1,000,000-1,000,500');
    expect(() => parseRegion(exampleRegion(seqs) as string, seqs)).not.toThrow();
  });

  it('falls back to the middle of a short sequence', () => {
    const seqs = [{ name: 'scaf_1', length: 2000 }];
    expect(exampleRegion(seqs)).toBe('scaf_1:750-1,250');
  });

  it('names a real sequence in the empty-input hint', () => {
    expect(() => parseRegion('', [{ name: 'SUPER_1', length: 3_000_000 }])).toThrow(
      /SUPER_1:1,000,000-1,000,500/,
    );
  });

  it('returns null with no sequences', () => {
    expect(exampleRegion([])).toBeNull();
  });
});
