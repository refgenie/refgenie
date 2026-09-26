/**
 * Flat lineage rows -> the nested taxonomy the radial layout consumes.
 *
 * One grouping pass per level, not a `.filter` per unique taxon per level, and
 * a closed row type instead of an index signature of `any`.
 */

import { childLevel, TAXONOMIC_LEVELS } from '../../types/taxonomy';
import type { TaxaRow, TaxonomicLevel, TreeNode } from '../../types/taxonomy';

function subtree(rows: TaxaRow[], level: TaxonomicLevel): TreeNode[] {
  const groups = new Map<string, TaxaRow[]>();
  for (const row of rows) {
    const name = row[level];
    if (!name) continue;
    const bucket = groups.get(name);
    if (bucket) bucket.push(row);
    else groups.set(name, [row]);
  }

  const next = childLevel(level);
  return [...groups.entries()].map(([name, members]) => ({
    name,
    taxonomicLevel: level,
    children: next ? subtree(members, next) : null,
  }));
}

/** The taxonomy hierarchy under a synthetic root. Pure; safe at module scope. */
export function buildTree(rows: TaxaRow[]): TreeNode {
  return { name: 'root', children: subtree(rows, TAXONOMIC_LEVELS[0]) };
}
