/** External-link derivation for genome metadata. */

/**
 * Taxonomy URI: the server-derived `taxon_uri` when present, otherwise the
 * identifiers.org form built from `taxon_id`.
 */
export function taxonUrl(
  taxonUri: string | null | undefined,
  taxonId: number | null | undefined,
): string | undefined {
  if (taxonUri) return taxonUri;
  if (taxonId === null || taxonId === undefined) return undefined;
  return `https://identifiers.org/taxonomy:${taxonId}`;
}

/**
 * NCBI datasets URL for a GCA_/GCF_ assembly accession. An FHR
 * `accessionID.url` wins when the sidecar supplies one.
 */
export function accessionUrl(
  accession: string | null | undefined,
  fhrAccessionUrl?: string | null,
): string | undefined {
  if (fhrAccessionUrl) return fhrAccessionUrl;
  if (!accession) return undefined;
  if (!/^GC[AF]_/.test(accession)) return undefined;
  return `https://www.ncbi.nlm.nih.gov/datasets/genome/${accession}/`;
}
