import { describe, expect, it } from 'vitest';
import { buildTreeForSpecies } from './taxa';
import type { TreeNode } from '../../types/taxonomy';

function leaves(node: TreeNode): string[] {
  if (!node.children || node.children.length === 0) return [node.name];
  return node.children.flatMap(leaves);
}

describe('buildTreeForSpecies', () => {
  it('keeps only the named species and reports the rest as missing', () => {
    const { tree, missing } = buildTreeForSpecies(['Homo sapiens', 'Not a species']);
    expect(leaves(tree).map((n) => n.toLowerCase())).toEqual(['homo sapiens']);
    expect(missing).toEqual(['Not a species']);
  });
});
