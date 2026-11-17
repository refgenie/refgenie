import { DownloadLink } from '../common/DownloadLink';
import { assetFileDownloadUrl } from '../../services/resources/serverInfo';
import { useApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';

/** Up to this asset size, files are renamed on download, which reads them into memory. */
const RENAME_MAX_BYTES = 50_000_000;

export interface AssetFileListProps {
  assetDigest: string;
  files: string[] | undefined;
  /** The asset's total size in bytes. */
  size?: number | null;
  /** The asset's genome, to save `<digest>.chrom.sizes` as `hg38.chrom.sizes`. */
  genome?: { digest: string; name: string };
}

/**
 * Download links for a file-mode asset's files, gated on
 * `capabilities.downloads`. Local mode uses `LocalAssetFileList` instead.
 */
export function AssetFileList({ assetDigest, files, size, genome }: AssetFileListProps) {
  const client = useApiClient();
  const canDownload = useCapability('downloads');
  const rename = !!genome && size != null && size <= RENAME_MAX_BYTES;
  const savedName = (file: string) =>
    rename && file.startsWith(genome.digest)
      ? genome.name + file.slice(genome.digest.length)
      : undefined;

  if (!files || files.length === 0) {
    return <p className="rg-muted text-sm">This asset is not staged for file-level serving.</p>;
  }

  return (
    <ul className="flex flex-col gap-1">
      {files.map((file) => (
        <li className="text-sm" key={file}>
          {canDownload ? (
            <DownloadLink
              href={assetFileDownloadUrl(client, assetDigest, file)}
              saveAs={savedName(file)}
            >
              <code className="rg-code rg-code--inline">{savedName(file) ?? file}</code>
            </DownloadLink>
          ) : (
            <code className="rg-code rg-code--inline">{file}</code>
          )}
        </li>
      ))}
    </ul>
  );
}
