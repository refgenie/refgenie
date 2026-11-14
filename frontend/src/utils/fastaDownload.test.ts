import { describe, expect, it } from 'vitest';
import {
  approxBytes,
  estimateFastaGzBytes,
  primaryChromosomes,
  refgetExportCommand,
  shellQuote,
} from './fastaDownload';

describe('primaryChromosomes', () => {
  it('drops unplaced, alternate, patch, decoy, EBV and HLA sequences', () => {
    const names = [
      'chr1',
      'chr2',
      'chrX',
      'chrY',
      'chrM',
      'chr1_KI270706v1_random',
      'chrUn_GL000220v1',
      'chr6_GL000250v2_alt',
      'chr1_KN196472v1_fix',
      'chrUn_JTFH01000001v1_decoy',
      'chrEBV',
      'HLA-A*01:01:01:01',
    ];
    expect(primaryChromosomes(names)).toEqual(['chr1', 'chr2', 'chrX', 'chrY', 'chrM']);
  });

  it('handles Ensembl and NCBI style names', () => {
    expect(primaryChromosomes(['1', '2', 'MT', 'KI270728.1', 'GL000009.2'])).toEqual([
      '1',
      '2',
      'MT',
    ]);
  });

  it('gives null when the pattern does not split the genome', () => {
    expect(primaryChromosomes(['SUPER_1', 'SUPER_2', 'scaffold_3'])).toBeNull();
    expect(primaryChromosomes(['chr1_random'])).toBeNull();
  });
});

describe('refgetExportCommand', () => {
  it('builds the whole-genome command', () => {
    expect(
      refgetExportCommand({
        digest: 'EiFob05a',
        storeUrl: 'https://s3.test/store/',
        output: 'hg38.fa.gz',
      }),
    ).toBe(
      'refget store export EiFob05a --remote https://s3.test/store/ -o hg38.fa.gz',
    );
  });

  it('adds -n per chosen sequence and quotes odd names', () => {
    expect(
      refgetExportCommand({
        digest: 'd',
        storeUrl: 'u',
        output: 'my genome.fa.gz',
        names: ['chr1', "it's*"],
      }),
    ).toBe(
      `refget store export d --remote u -n chr1 -n 'it'\\''s*' -o 'my genome.fa.gz'`,
    );
  });

  it('leaves plain words unquoted', () => {
    expect(shellQuote('chr1_KI270706v1_random')).toBe('chr1_KI270706v1_random');
  });
});

describe('estimateFastaGzBytes and approxBytes', () => {
  it('puts hg38 at about 0.9 GB', () => {
    expect(approxBytes(estimateFastaGzBytes(3_099_922_541))).toBe('0.9 GB');
  });
  it('rounds to two figures, moving up a unit from half of it', () => {
    expect(approxBytes(300)).toBe('300 B');
    expect(approxBytes(600)).toBe('0.6 KB');
    expect(approxBytes(12_345)).toBe('12 KB');
    expect(approxBytes(300 * 1024 ** 2)).toBe('300 MB');
    expect(approxBytes(640 * 1024 ** 2)).toBe('0.6 GB');
    expect(approxBytes(2.5 * 1024 ** 3)).toBe('2.5 GB');
  });
});
