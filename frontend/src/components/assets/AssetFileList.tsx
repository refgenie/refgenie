import { assetFileDownloadUrl } from '../../services/resources/serverInfo';
import { useApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';

export interface AssetFileListProps {
  assetDigest: string;
  files: string[] | undefined;
}

/**
 * The individually-downloadable files of a file-mode asset. Download links are
 * real anchors built by `ApiClient.url()`, gated on `capabilities.downloads`
 * (the local dash has no file-download route).
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
            <a className="rg-link" href={assetFileDownloadUrl(client, assetDigest, file)}>
              <code className="rg-code rg-code--inline">{file}</code>
            </a>
          ) : (
            <code className="rg-code rg-code--inline">{file}</code>
          )}
        </li>
      ))}
    </ul>
  );
}
