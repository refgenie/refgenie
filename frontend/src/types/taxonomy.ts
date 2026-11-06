/** Taxonomy types for the tree explorer, built from `src/assets/taxa.json` (not the API). */

export const TAXONOMIC_LEVELS = [
  'domain',
  'kingdom',
  'phylum',
  'class',
  'order',
  'family',
  'genus',
  'species',
] as const;

export type TaxonomicLevel = (typeof TAXONOMIC_LEVELS)[number];

/** One row of `taxa.json`: a fully-qualified lineage. */
export type TaxaRow = Record<TaxonomicLevel, string>;

export interface TreeNode {
  name: string;
  /** Absent only on the synthetic root. */
  taxonomicLevel?: TaxonomicLevel;
  children?: TreeNode[] | null;
}

/** The rank below `level`, or `null` at the leaves. */
export function childLevel(level: TaxonomicLevel): TaxonomicLevel | null {
  const i = TAXONOMIC_LEVELS.indexOf(level);
  return i < 0 || i === TAXONOMIC_LEVELS.length - 1 ? null : TAXONOMIC_LEVELS[i + 1];
}
