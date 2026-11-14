/** External-link derivation for genome metadata. */

/**
 * `url` if it is an http(s) link, else undefined. Relative links resolve
 * against the page, so they pass. Use it on every href or src that comes from
 * API, metadata or config data: a `javascript:` URL would run on this site.
 */
export function safeHttpUrl(url: string | null | undefined): string | undefined {
  if (typeof url !== 'string' || url.trim() === '') return undefined;
  const base = typeof document !== 'undefined' ? document.baseURI : undefined;
  let parsed: URL;
  try {
    parsed = new URL(url, base);
  } catch {
    return undefined;
  }
  return parsed.protocol === 'http:' || parsed.protocol === 'https:' ? url : undefined;
}

/**
 * Taxonomy URI: the server-derived `taxon_uri` when present, otherwise the
 * identifiers.org form built from `taxon_id`.
 */
export function taxonUrl(
  taxonUri: string | null | undefined,
  taxonId: number | null | undefined,
): string | undefined {
  const safe = safeHttpUrl(taxonUri);
  if (safe) return safe;
  if (taxonId === null || taxonId === undefined) return undefined;
  return `https://identifiers.org/taxonomy:${taxonId}`;
}

/**
 * NCBI datasets URL for a GCA_/GCF_ assembly accession. An FHR
 * `accessionID.url` wins when the sidecar supplies an http(s) one.
 */
export function accessionUrl(
  accession: string | null | undefined,
  fhrAccessionUrl?: string | null,
): string | undefined {
  const safe = safeHttpUrl(fhrAccessionUrl);
  if (safe) return safe;
  if (!accession) return undefined;
  if (!/^GC[AF]_/.test(accession)) return undefined;
  return `https://www.ncbi.nlm.nih.gov/datasets/genome/${accession}/`;
}
