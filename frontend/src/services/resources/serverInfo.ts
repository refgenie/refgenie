/**
 * Server-mode-only reads and download URLs. Gate calls on
 * `capabilities.archives` / `capabilities.downloads`, not on the mode.
 */

import { listQuery } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { ArchiveRecord, DatabaseSummaryResponse } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

export const getSummary = (c: ApiClient, init?: RequestInitLite) =>
  c.get<DatabaseSummaryResponse>('/summary', undefined, init);

export interface ListArchivesParams extends ListParams {
  genome_digest?: string;
  asset_digest?: string;
}

export const listArchives = (
  c: ApiClient,
  p: ListArchivesParams = {},
  init?: RequestInitLite,
) =>
  c.get<PageResponse<ArchiveRecord>>(
    '/archives',
    {
      genome_digest: p.genome_digest,
      asset_digest: p.asset_digest,
      ...listQuery(p),
    },
    init,
  );

/** `<a href>` target for the archive tarball of an asset. */
export const archiveDownloadUrl = (c: ApiClient, assetDigest: string): string =>
  c.url(`/archives/${encodeURIComponent(assetDigest)}/download`);

/** `<a href>` target for one file in an asset; `/` separators stay unencoded. */
export const assetFileDownloadUrl = (
  c: ApiClient,
  assetDigest: string,
  filePath: string,
): string => {
  const encoded = filePath
    .split('/')
    .map((segment) => encodeURIComponent(segment))
    .join('/');
  return c.url(`/assets/${encodeURIComponent(assetDigest)}/files/${encoded}`);
};
