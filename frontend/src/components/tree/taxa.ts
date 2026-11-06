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
const CANONICAL = new Map(TAXA.map((row) => [normalizeSpecies(row.species), row.species]));

/** The `taxa.json` spelling of `name`, or null if the tree does not have it. */
export function resolveSpecies(name: string): string | null {
  return CANONICAL.get(normalizeSpecies(name)) ?? null;
}
