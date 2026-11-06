import { describe, expect, it } from 'vitest';
import { buildTree } from './buildTree';
import { TAXONOMIC_LEVELS } from '../../types/taxonomy';
import type { TaxaRow, TreeNode } from '../../types/taxonomy';

const row = (over: Partial<TaxaRow>): TaxaRow => ({
  domain: 'Eukaryota',
  kingdom: 'Metazoa',
  phylum: 'Chordata',
  class: 'Mammalia',
  order: 'Primates',
  family: 'Hominidae',
  genus: 'Homo',
  species: 'Homo sapiens',
  ...over,
});

/** Walk `depth` levels down the first child at each step. */
function descend(node: TreeNode | undefined, depth: number): TreeNode | undefined {
  let current = node;
  for (let i = 0; i < depth; i += 1) current = current?.children?.[0];
  return current;
}

describe('buildTree', () => {
  it('nests one lineage eight levels deep', () => {
    let node = buildTree([row({})]).children?.[0];
    for (const level of TAXONOMIC_LEVELS) {
      expect(node?.taxonomicLevel).toBe(level);
      node = node?.children?.[0];
    }
    expect(node).toBeUndefined();
  });

  it('merges rows that share a prefix and splits where they diverge', () => {
    const root = buildTree([row({}), row({ genus: 'Pan', species: 'Pan troglodytes' })]);
    // root -> domain -> kingdom -> phylum -> class -> order -> family
    const family = descend(root, 6);
    expect(family?.name).toBe('Hominidae');
    expect(family?.children?.map((c) => c.name)).toEqual(['Homo', 'Pan']);
  });

  it('leaves species children null', () => {
    const species = descend(buildTree([row({})]), 8);
    expect(species?.name).toBe('Homo sapiens');
    expect(species?.children).toBeNull();
  });
});
