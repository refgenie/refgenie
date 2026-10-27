import { useState } from 'react';
import { DownloadLink } from '../common/DownloadLink';
import { FileSize } from '../common/FileSize';
import { ErrorState, LoadingState } from '../common/states';
import { LocalFilePreviewModal } from './LocalFilePreviewModal';
import { useLocalAssetFiles } from '../../hooks/queries/useLocalFiles';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { localFileDownloadUrl } from '../../services/resources/localFiles';

export interface LocalAssetFileListProps {
  assetDigest: string;
}

/**
 * Every file in a local asset's directory, read from disk through `/v1/files`.
 * Local mode only: the caller gates it on `capabilities.local_files`.
 *
 * Each file name downloads the file. Text files also get a "View" button that
 * opens a capped preview; binary and compressed files are download-only.
 * Directories are listed as their own rows, so nested index files read as a
 * tree of paths.
 */
export function LocalAssetFileList({ assetDigest }: LocalAssetFileListProps) {
  const client = useLocalApiClient();
  const files = useLocalAssetFiles(assetDigest);
  const [previewPath, setPreviewPath] = useState<string | null>(null);

  if (files.isPending) return <LoadingState label="Loading files" />;
  if (files.error) {
    return <ErrorState error={files.error} subject={assetDigest} />;
  }

  const entries = files.data?.entries ?? [];
  if (entries.length === 0) {
    return <p className="rg-muted text-sm">This asset has no files on disk.</p>;
  }

  return (
    <>
      <ul className="flex flex-col gap-1">
        {entries.map((entry) => (
          <li className="flex items-center gap-2 flex-wrap text-sm" key={entry.path}>
            {entry.type === 'dir' ? (
              <code className="rg-code rg-code--inline rg-muted">{entry.path}/</code>
            ) : (
              <>
                <DownloadLink href={localFileDownloadUrl(client, assetDigest, entry.path)}>
                  <code className="rg-code rg-code--inline">{entry.path}</code>
                </DownloadLink>
                <span className="rg-muted text-xs">
                  <FileSize bytes={entry.size} />
                </span>
                {entry.previewable && (
                  <button
                    type="button"
                    className="rg-btn rg-btn--sm"
                    onClick={() => setPreviewPath(entry.path)}
                    aria-label={`View ${entry.path}`}
                  >
                    View
                  </button>
                )}
              </>
            )}
          </li>
        ))}
      </ul>
      {files.data?.truncated && (
        <p className="rg-muted text-sm mt-2">
          Only the first {entries.length} entries are listed.
        </p>
      )}
      <LocalFilePreviewModal
        assetDigest={assetDigest}
        path={previewPath}
        onClose={() => setPreviewPath(null)}
      />
    </>
  );
}
