import { describe, expect, it } from 'vitest';
import { buildGenomeIndex, genomeLabel, genomeMatches, preferredAlias } from './genomes';
import type { GenomeResponse } from '../types/api';

const DIGEST = '5cb5b9d2d1e2b5e5e2f0a3c4d5e6f708';

function genome(overrides: Partial<GenomeResponse> = {}): GenomeResponse {
  return {
    digest: DIGEST,
    aliases: ['hg38'],
    description: null,
    asset_count: 0,
    species_name: null,
    common_name: null,
    taxon_id: null,
    assembly_source: null,
    assembly_accession: null,
    ...overrides,
  };
}

/** The live alias array of the canonical human genome on api.refgenie.org. */
const HUMAN = [
  'hg38-refgenie',
  'hg38-noALT-noHLA-noDecoy-broad',
  'hg38',
  'GRCh38.p14-fasta-no-alt-analysis',
  'hs38',
  'GRCh38-igenomes-ncbi',
];

describe('preferredAlias', () => {
  it('echoes the term the user searched for rather than aliases[0]', () => {
    expect(preferredAlias(HUMAN, 'hg38')).toBe('hg38');
  });

  it('falls back to a substring match when nothing matches exactly', () => {
    expect(preferredAlias(HUMAN, 'igenomes')).toBe('GRCh38-igenomes-ncbi');
  });

  it('is case-insensitive and ignores surrounding whitespace', () => {
    expect(preferredAlias(HUMAN, '  HG38 ')).toBe('hg38');
  });

  it('picks the shortest alias when there is no query, ties broken alphabetically', () => {
    expect(preferredAlias(HUMAN)).toBe('hg38');
  });

  it('breaks a length tie alphabetically, so the choice is stable', () => {
    expect(preferredAlias(['mm10', 'hg38'])).toBe('hg38');
  });

  it('returns undefined for an empty array', () => {
    expect(preferredAlias([])).toBeUndefined();
    expect(preferredAlias(undefined)).toBeUndefined();
  });

  it('does not mutate the array it was given', () => {
    const aliases = [...HUMAN];
    preferredAlias(aliases);
    expect(aliases).toEqual(HUMAN);
  });
});

describe('genomeLabel', () => {
  it('returns the preferred alias', () => {
    const index = buildGenomeIndex([genome({ aliases: ['hg38', 'GRCh38'] })]);
    expect(genomeLabel(index, DIGEST)).toBe('hg38');
  });

  it('leads with the alias matching the query', () => {
    const index = buildGenomeIndex([genome({ aliases: HUMAN })]);
    expect(genomeLabel(index, DIGEST, 'GRCh38.p14')).toBe(
      'GRCh38.p14-fasta-no-alt-analysis',
    );
  });

  it('falls back to a truncated digest for an unknown digest', () => {
    const index = buildGenomeIndex([genome()]);
    expect(genomeLabel(index, 'a1b2c3d4e5f60718293a4b5c6d7e8f90')).toBe('a1b2c3d4e5f6…');
  });

  it('falls back to a truncated digest when the genome carries no alias', () => {
    const index = buildGenomeIndex([genome({ aliases: [] })]);
    expect(genomeLabel(index, DIGEST)).toBe('5cb5b9d2d1e2…');
  });

  it('returns NA for a null digest', () => {
    expect(genomeLabel(buildGenomeIndex([]), null)).toBe('NA');
  });

  it('does not throw while the index is still loading', () => {
    expect(genomeLabel(undefined, DIGEST)).toBe('5cb5b9d2d1e2…');
  });
});

describe('genomeMatches', () => {
  const hg38 = genome({
    aliases: ['hg38', 'GRCh38'],
    description: 'Human GRCh38 reference assembly',
    species_name: 'Homo sapiens',
    assembly_accession: 'GCA_000001405.15',
  });

  it('matches an empty query', () => {
    expect(genomeMatches(hg38, '', [], 'contains')).toBe(true);
    expect(genomeMatches(hg38, '   ', [], 'contains')).toBe(true);
  });

  it('searches every field when none is selected', () => {
    expect(genomeMatches(hg38, 'GRCh38', [], 'contains')).toBe(true);
    expect(genomeMatches(hg38, 'sapiens', [], 'contains')).toBe(true);
    expect(genomeMatches(hg38, 'GCA_0000', [], 'contains')).toBe(true);
    expect(genomeMatches(hg38, 'mm10', [], 'contains')).toBe(false);
  });

  it('honours the selected field, so the field selector is not decorative', () => {
    expect(genomeMatches(hg38, 'sapiens', ['species_name'], 'contains')).toBe(true);
    expect(genomeMatches(hg38, 'sapiens', ['aliases'], 'contains')).toBe(false);
  });

  it('honours the operator, matching the server: eq is case-sensitive', () => {
    expect(genomeMatches(hg38, 'hg38', ['aliases'], 'eq')).toBe(true);
    expect(genomeMatches(hg38, 'HG38', ['aliases'], 'eq')).toBe(false);
    expect(genomeMatches(hg38, 'HG', ['aliases'], 'starts_with')).toBe(true);
    expect(genomeMatches(hg38, '38', ['aliases'], 'ends_with')).toBe(true);
    expect(genomeMatches(hg38, 'hg', ['aliases'], 'ends_with')).toBe(false);
  });
});
