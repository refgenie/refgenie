/**
 * The local digest index: every genome and asset digest present on the
 * connected local refgenie, fetched once and cached. Presence badges are an
 * exact intersection of these digests with the remote listing's digests —
 * digests are content-derived and identical across instances, which is the
 * whole reason this feature is worth building (no fuzzy name matching).
 */

import { listAssets } from '../resources/assets';
import { listGenomes } from '../resources/genomes';
import type { ApiClient } from '../http';

/** The server cap (`refgenie/const.py::MAX_PAGE_SIZE`); over it is a 422. */
const PAGE_SIZE = 1000;

/** A local instance is small; this only stops a pathological one looping. */
const MAX_PAGES = 20;

export interface BridgeDigestIndex {
  genomeDigests: Set<string>;
  assetDigests: Set<string>;
}

async function collect(
  fetchPage: (offset: number) => Promise<{
    items: Array<{ digest: string | null }>;
    pagination: { total: number };
  }>,
): Promise<Set<string>> {
  const digests = new Set<string>();
  let offset = 0;
  for (let page = 0; page < MAX_PAGES; page++) {
    const data = await fetchPage(offset);
    for (const item of data.items) {
      // AssetResponse.digest is genuinely nullable on the wire.
      if (typeof item.digest === 'string') digests.add(item.digest);
    }
    offset += data.items.length;
    if (data.items.length < PAGE_SIZE || offset >= data.pagination.total) break;
  }
  return digests;
}

export async function fetchBridgeDigestIndex(
  read: ApiClient,
): Promise<BridgeDigestIndex> {
  const [genomeDigests, assetDigests] = await Promise.all([
    collect((offset) => listGenomes(read, { offset, limit: PAGE_SIZE })),
    collect((offset) => listAssets(read, { offset, limit: PAGE_SIZE })),
  ]);
  return { genomeDigests, assetDigests };
}
