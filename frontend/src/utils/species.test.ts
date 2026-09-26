import { describe, expect, it } from 'vitest';
import {
  buildSpeciesIndex,
  compareSpecies,
  normalizeSpecies,
  speciesGenomesHref,
  speciesMatches,
} from './species';
import type { SpeciesRow } from './species';
import type { GenomeResponse } from '../types/api';

let counter = 0;

function genome(overrides: Partial<GenomeResponse> = {}): GenomeResponse {
  counter += 1;
  return {
    digest: `digest-${counter}`,
    aliases: [],
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

const HUMAN = { species_name: 'Homo sapiens', common_name: 'human', taxon_id: 9606 };

function row(overrides: Partial<SpeciesRow> = {}): SpeciesRow {
  return {
    key: 'homo sapiens',
    speciesName: 'Homo sapiens',
    commonNames: ['human'],
    taxonIds: [9606],
    genomes: 1,
    assets: 0,
    ...overrides,
  };
}

describe('buildSpeciesIndex', () => {
  it('groups genomes of one species and sums their counts', () => {
    const index = buildSpeciesIndex([
      genome({ ...HUMAN, asset_count: 3 }),
      genome({ ...HUMAN, asset_count: 1 }),
    ]);
    expect(index.size).toBe(1);
    const human = index.get('homo sapiens')!;
    expect(human.genomes).toBe(2);
    expect(human.assets).toBe(4);
  });

  it('folds the taxonomy spelling into the API spelling, and renders the majority', () => {
    // `taxa.json` says "Homo Sapiens"; the API says "Homo sapiens". They
    // must fold into one row.
    const index = buildSpeciesIndex([
      genome({ ...HUMAN, species_name: 'Homo sapiens' }),
      genome({ ...HUMAN, species_name: 'Homo sapiens' }),
      genome({ ...HUMAN, species_name: 'Homo Sapiens' }),
    ]);
    expect(index.size).toBe(1);
    expect(index.get('homo sapiens')!.speciesName).toBe('Homo sapiens');
    expect(index.get('homo sapiens')!.genomes).toBe(3);
  });

  it('collects distinct common names, sorted and deduped', () => {
    const index = buildSpeciesIndex([
      genome({ ...HUMAN, common_name: 'human' }),
      genome({ ...HUMAN, common_name: 'human' }),
      genome({ ...HUMAN, common_name: 'Human being' }),
    ]);
    expect(index.get('homo sapiens')!.commonNames).toEqual(['human', 'Human being']);
  });

  it('collects distinct taxon IDs and skips nulls', () => {
    const index = buildSpeciesIndex([
      genome({ ...HUMAN, taxon_id: 9606 }),
      genome({ ...HUMAN, taxon_id: null }),
      genome({ ...HUMAN, taxon_id: 63221 }),
    ]);
    expect(index.get('homo sapiens')!.taxonIds).toEqual([9606, 63221]);
  });

  it('puts genomes with no species into the unspecified row', () => {
    const index = buildSpeciesIndex([genome({ asset_count: 2 }), genome({ ...HUMAN })]);
    const unspecified = index.get('')!;
    expect(unspecified.genomes).toBe(1);
    expect(unspecified.speciesName).toBe('');
    expect(speciesGenomesHref(unspecified)).toBeNull();
  });
});

describe('normalizeSpecies', () => {
  it('case-folds and collapses whitespace', () => {
    expect(normalizeSpecies('  Homo   Sapiens ')).toBe('homo sapiens');
    expect(normalizeSpecies(null)).toBe('');
  });
});

describe('speciesMatches', () => {
  const human = row();

  it('matches on the common name', () => {
    expect(speciesMatches(human, 'human')).toBe(true);
  });

  it('matches on the scientific name', () => {
    expect(speciesMatches(human, 'homo')).toBe(true);
  });

  it('matches on a bare taxon ID', () => {
    expect(speciesMatches(human, '9606')).toBe(true);
  });

  it('misses on an unrelated term', () => {
    expect(speciesMatches(human, 'mouse')).toBe(false);
  });

  it('matches everything on an empty query', () => {
    expect(speciesMatches(human, '   ')).toBe(true);
  });
});

describe('speciesGenomesHref', () => {
  it('searches the species column with a case-insensitive operator', () => {
    expect(speciesGenomesHref(row({ assets: 4 }))).toBe(
      '/genomes?q=Homo+sapiens&fields=species_name&op=contains',
    );
  });

  it('unticks the assets filter only when the species has nothing built', () => {
    expect(speciesGenomesHref(row({ assets: 0 }))).toBe(
      '/genomes?q=Homo+sapiens&fields=species_name&op=contains&assets=all',
    );
  });
});

describe('compareSpecies', () => {
  it('sorts the unspecified row last no matter how many genomes it holds', () => {
    const unspecified = row({ key: '', speciesName: '', genomes: 99, assets: 99 });
    const human = row({ assets: 1 });
    expect([unspecified, human].sort(compareSpecies)[0]).toBe(human);
    expect([human, unspecified].sort(compareSpecies)[1]).toBe(unspecified);
  });

  it('sorts by assets, then genomes, then name', () => {
    const many = row({ key: 'a', speciesName: 'A', assets: 10 });
    const few = row({ key: 'b', speciesName: 'B', assets: 1 });
    expect([few, many].sort(compareSpecies)).toEqual([many, few]);
  });
});
