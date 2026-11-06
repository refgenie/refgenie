/**
 * Local file access (local mode, `capabilities.local_files`) via the `/v1`
 * client. The server refuses these from other origins, including the bridge.
 */

import type { ApiClient, RequestInitLite } from '../http';
import type { LocalAssetFiles, LocalFilePreview } from '../../types/api';

/** Encode each path segment but keep the `/` separators (`{file_path:path}`). */
function encodePath(filePath: string): string {
  return filePath
    .split('/')
    .map((segment) => encodeURIComponent(segment))
    .join('/');
}

const assetBase = (assetDigest: string) => `/files/assets/${encodeURIComponent(assetDigest)}`;

export const listLocalAssetFiles = (c: ApiClient, assetDigest: string, init?: RequestInitLite) =>
  c.get<LocalAssetFiles>(assetBase(assetDigest), undefined, init);

export const getLocalFilePreview = (
  c: ApiClient,
  assetDigest: string,
  filePath: string,
  init?: RequestInitLite,
) => c.get<LocalFilePreview>(`${assetBase(assetDigest)}/preview/${encodePath(filePath)}`, undefined, init);

/** `<a href>` target that downloads one file (`Content-Disposition: attachment`). */
export const localFileDownloadUrl = (c: ApiClient, assetDigest: string, filePath: string): string =>
  c.url(`${assetBase(assetDigest)}/raw/${encodePath(filePath)}`);
