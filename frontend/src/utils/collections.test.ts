import { describe, expect, it } from 'vitest';
import {
  collectionStats,
  expandGenomeLink,
  filterCollectionGenomes,
} from './collections';
import { vgpGenomes } from '../test/collectionFixtures';

describe('collectionStats', () => {
  it('counts genomes, species, assets and haplotypes', () => {
    const stats = collectionStats(vgpGenomes);
    expect(stats.genomes).toBe(3);
    expect(stats.species).toBe(2);
    expect(stats.withAssets).toBe(1);
    expect(stats.assetGroups).toBe(2);
    expect(stats.haplotypes).toMatchObject({
      maternal: 1,
      paternal: 1,
      primary: 1,
      hap1: 0,
    });
  });
});

describe('expandGenomeLink', () => {
  const genome = vgpGenomes[0];

  it('fills every placeholder', () => {
    expect(
      expandGenomeLink(
        'https://x/{assembly_accession}/{species_underscored}/{species_name}/{digest}',
        genome,
      ),
    ).toBe(
      `https://x/GCF_902713425.1/Acipenser_ruthenus/Acipenser%20ruthenus/${genome.digest}`,
    );
  });

  it('returns null when a used placeholder is missing', () => {
    expect(
      expandGenomeLink('https://x/{assembly_accession}', {
        ...genome,
        assembly_accession: null,
      }),
    ).toBeNull();
    expect(expandGenomeLink('https://x/{nope}', genome)).toBeNull();
  });
});

describe('filterCollectionGenomes', () => {
  it('matches accession and alias, case-insensitively', () => {
    expect(filterCollectionGenomes(vgpGenomes, { q: 'gca_905165185' })).toEqual([
      vgpGenomes[1],
    ]);
    expect(filterCollectionGenomes(vgpGenomes, { q: 'FLATCHA' })).toEqual([
      vgpGenomes[2],
    ]);
  });

  it('respects the haplotype and assets filters', () => {
    expect(filterCollectionGenomes(vgpGenomes, { haplotype: 'paternal' })).toEqual([
      vgpGenomes[1],
    ]);
    expect(filterCollectionGenomes(vgpGenomes, { onlyWithAssets: true })).toEqual([
      vgpGenomes[0],
    ]);
    expect(
      filterCollectionGenomes(vgpGenomes, { q: 'acipenser', haplotype: 'primary' }),
    ).toEqual([vgpGenomes[2]]);
  });
});
