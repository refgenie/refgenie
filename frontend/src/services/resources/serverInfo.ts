/**
 * Server-mode-only reads (`refgenie/server/routers/version4.py`) plus the
 * download-URL builders.
 *
 * Gate calls on `capabilities.archives` / `capabilities.downloads`, never on
 * the mode string.
 */

import type { ApiClient, RequestInitLite } from '../http';
import type {
  ArchiveRecord,
  DatabaseSummaryResponse,
  SpeciesSummaryResponse,
} from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

export const getSummary = (c: ApiClient, init?: RequestInitLite) =>
  c.get<DatabaseSummaryResponse>('/summary', undefined, init);

export const getSpeciesSummary = (c: ApiClient, init?: RequestInitLite) =>
  c.get<SpeciesSummaryResponse>('/species/summary', undefined, init);

export interface ListArchivesParams extends ListParams {
  genome_digest?: string;
  asset_digest?: string;
}

export const listArchives = (
  c: ApiClient,
  p: ListArchivesParams = {},
  init?: RequestInitLite,
) =>
  c.get<Paginated<ArchiveRecord>>(
    '/archives',
    {
      genome_digest: p.genome_digest,
      asset_digest: p.asset_digest,
      offset: p.offset,
      limit: p.limit,
    },
    init,
  );

/** `<a href>` target for the archive tarball of an asset. */
export const archiveDownloadUrl = (c: ApiClient, assetDigest: string): string =>
  c.url(`/archives/${encodeURIComponent(assetDigest)}/download`);

/**
 * `<a href>` target for one file inside an asset.
 *
 * The route is `{file_path:path}`, so each segment is encoded but the `/`
 * separators are left intact.
 */
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
