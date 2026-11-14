/**
 * Reads the bundled `taxa.json` (~130 KB). Import it only from the lazy `/tree`
 * chunk; `/species` gets its species from the API.
 */

import taxaRaw from '../../assets/taxa.json';
import { buildTree } from './buildTree';
import { normalizeSpecies } from '../../utils/species';
import type { TaxaRow, TreeNode } from '../../types/taxonomy';

const TAXA = taxaRaw as TaxaRow[];

/** Pure and deterministic, so it is built once at module load, not per render. */
export const treeData: TreeNode = buildTree(TAXA);

/**
 * Normalized species name -> the spelling `taxa.json` uses.
 *
 * `taxa.json` says "Homo Sapiens"; the API says "Homo sapiens". `Tree` selects
 * a tip with `d.data.name === selected`, so without this fold an inbound
 * `/tree?species=Homo sapiens` link highlights nothing.
 */
const CANONICAL = new Map(
  TAXA.map((row) => [normalizeSpecies(row.species), row.species]),
);

/** The `taxa.json` spelling of `name`, or null if the tree does not have it. */
export function resolveSpecies(name: string): string | null {
  return CANONICAL.get(normalizeSpecies(name)) ?? null;
}

/**
 * The taxonomy cut down to these species (any spelling), plus the names it
 * does not have. Kept here so `taxa.json` stays in the lazy chunk.
 */
export function buildTreeForSpecies(names: Iterable<string>): {
  tree: TreeNode;
  missing: string[];
} {
  const wanted = new Map<string, string>();
  for (const name of names) {
    const key = normalizeSpecies(name);
    if (key && !wanted.has(key)) wanted.set(key, name);
  }
  const rows = TAXA.filter((row) => wanted.has(normalizeSpecies(row.species)));
  const found = new Set(rows.map((row) => normalizeSpecies(row.species)));
  const missing = [...wanted.entries()]
    .filter(([key]) => !found.has(key))
    .map(([, name]) => name);
  return { tree: buildTree(rows), missing };
}
