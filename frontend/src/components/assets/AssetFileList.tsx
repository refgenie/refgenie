import { DownloadLink } from '../common/DownloadLink';
import { assetFileDownloadUrl } from '../../services/resources/serverInfo';
import { useApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';

export interface AssetFileListProps {
  assetDigest: string;
  files: string[] | undefined;
}

/**
 * Download links for a file-mode asset's files, gated on
 * `capabilities.downloads`. Local mode uses `LocalAssetFileList` instead.
 */
export function AssetFileList({ assetDigest, files }: AssetFileListProps) {
  const client = useApiClient();
  const canDownload = useCapability('downloads');

  if (!files || files.length === 0) {
    return <p className="rg-muted text-sm">This asset is not staged for file-level serving.</p>;
  }

  return (
    <ul className="flex flex-col gap-1">
      {files.map((file) => (
        <li className="text-sm" key={file}>
          {canDownload ? (
            <DownloadLink href={assetFileDownloadUrl(client, assetDigest, file)}>
              <code className="rg-code rg-code--inline">{file}</code>
            </DownloadLink>
          ) : (
            <code className="rg-code rg-code--inline">{file}</code>
          )}
        </li>
      ))}
    </ul>
  );
}
