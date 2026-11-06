/**
 * Every genome and asset digest on the connected local refgenie. Digests are
 * content-derived, so "already local" badges are an exact match.
 */

import { listAssets } from '../resources/assets';
import { listGenomes } from '../resources/genomes';
import { walkPages } from '../../utils/pagination';
import type { ApiClient } from '../http';
import type { PageResponse } from '../../types/pagination';

/** The server cap (`refgenie/const.py::MAX_PAGE_SIZE`); over it is a 400. */
const PAGE_SIZE = 1000;

/** A local instance is small; this only stops a pathological one looping. */
const MAX_ITEMS = 20_000;

export interface BridgeDigestIndex {
  genomeDigests: Set<string>;
  assetDigests: Set<string>;
}

async function collect(
  fetchPage: (
    page: number,
    pageSize: number,
  ) => Promise<PageResponse<{ digest: string | null }>>,
): Promise<Set<string>> {
  const items = await walkPages(fetchPage, {
    pageSize: PAGE_SIZE,
    maxItems: MAX_ITEMS,
  });
  // AssetResponse.digest can be null.
  return new Set(
    items.map((item) => item.digest).filter((d): d is string => typeof d === 'string'),
  );
}

export async function fetchBridgeDigestIndex(
  read: ApiClient,
): Promise<BridgeDigestIndex> {
  const [genomeDigests, assetDigests] = await Promise.all([
    collect((page, pageSize) => listGenomes(read, { page, pageSize })),
    collect((page, pageSize) => listAssets(read, { page, pageSize })),
  ]);
  return { genomeDigests, assetDigests };
}
