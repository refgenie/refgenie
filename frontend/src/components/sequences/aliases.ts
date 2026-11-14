// Ported from repos/refget/frontend/src/services/aliases.ts @ c0d3e09; keep the two in step.

import type { AliasEntry } from './types';

/**
 * Display preference when a sequence or collection has aliases in several
 * namespaces. Purely a UI choice.
 */
export const COLLECTION_NAMESPACE_PRIORITY = [
  'genome_assembly',
  'name',
  'accession',
  'refseq',
  'insdc',
];

/**
 * Choose one alias to show, preferring COLLECTION_NAMESPACE_PRIORITY and falling
 * back to whatever the store carries. `all` is a deduped "namespace: alias"
 * list, suitable for a tooltip.
 */
export function preferredAlias(entries: AliasEntry[] | undefined | null): {
  primary: string | null;
  all: string[];
} {
  if (!entries || entries.length === 0) return { primary: null, all: [] };
  const all = [...new Set(entries.map((e) => `${e.namespace}: ${e.alias}`))];
  for (const namespace of COLLECTION_NAMESPACE_PRIORITY) {
    const hit = entries.find((e) => e.namespace === namespace);
    if (hit) return { primary: hit.alias, all };
  }
  return { primary: entries[0].alias, all };
}
