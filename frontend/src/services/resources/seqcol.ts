/**
 * The server's GA4GH seqcol APIs, one per store (`seqcol.services` in
 * service-info). Used to link a genome's digests to the API and to find other
 * genomes that share one of its attributes.
 */

import { ApiClient } from '../http';
import type { RequestInitLite } from '../http';
import { instanceUrl } from '../instanceLinks';
import type { SeqcolServiceInfo } from '../../types/ui';

/** The seqcol service that reads this store, else the default (first) one. */
export function serviceForStore(
  services: SeqcolServiceInfo[],
  storeUrl: string | null | undefined,
): SeqcolServiceInfo | undefined {
  const norm = (u: string) => u.replace(/\/+$/, '');
  const match = storeUrl
    ? services.find(
        (s) => s.refget_store_url && norm(s.refget_store_url) === norm(storeUrl),
      )
    : undefined;
  return match ?? services[0];
}

/** Absolute-or-root-relative base of a seqcol service, e.g. `https://api…/seqcol/vgp`. */
export function seqcolBase(apiBaseUrl: string, service: SeqcolServiceInfo): string {
  return /^https?:\/\//.test(service.url)
    ? service.url
    : instanceUrl(apiBaseUrl, service.url);
}

/** `GET {base}/collection/{digest}?level=1`: the collection's level-1 digests. */
export function collectionUrl(base: string, digest: string): string {
  return `${base}/collection/${encodeURIComponent(digest)}?level=1`;
}

/** `GET {base}/attribute/collection/{attribute}/{digest}`: the attribute's value. */
export function attributeUrl(base: string, attribute: string, digest: string): string {
  return `${base}/attribute/collection/${encodeURIComponent(attribute)}/${encodeURIComponent(digest)}`;
}

interface ListResponse {
  results: string[];
}

/**
 * Digests of every collection, in any of the server's stores, whose
 * `attribute` has this digest. Stores that fail to answer are skipped.
 */
export async function findCollectionsSharing(
  apiBaseUrl: string,
  services: SeqcolServiceInfo[],
  attribute: string,
  digest: string,
  init?: RequestInitLite,
): Promise<string[]> {
  const lists = await Promise.allSettled(
    services.map((service) =>
      new ApiClient({ baseUrl: seqcolBase(apiBaseUrl, service) }).get<ListResponse>(
        '/list/collection',
        { [attribute]: digest, page_size: 100 },
        init,
      ),
    ),
  );
  const found = new Set<string>();
  for (const list of lists) {
    if (list.status === 'fulfilled') list.value.results.forEach((d) => found.add(d));
  }
  return [...found];
}
