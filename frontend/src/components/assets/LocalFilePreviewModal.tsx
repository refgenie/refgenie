import { BaseModal } from '../common/BaseModal';
import { DownloadLink } from '../common/DownloadLink';
import { ErrorState } from '../common/states';
import { useLocalFilePreview } from '../../hooks/queries/useLocalFiles';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { localFileDownloadUrl } from '../../services/resources/localFiles';
import { formatBytes } from '../../utils/format';

export interface LocalFilePreviewModalProps {
  assetDigest: string;
  /** The file to preview; `null` closes the modal. */
  path: string | null;
  onClose: () => void;
}

/** The start of one text file from a local asset (server caps it at 64 KB). */
export function LocalFilePreviewModal({ assetDigest, path, onClose }: LocalFilePreviewModalProps) {
  const client = useLocalApiClient();
  const preview = useLocalFilePreview(assetDigest, path);
  const data = preview.data;

  return (
    <BaseModal isOpen={path !== null} onClose={onClose} title={path ?? ''} size="lg">
      <div className="flex flex-col gap-4">
        {preview.isPending ? (
          <p className="rg-muted text-sm">Loading preview…</p>
        ) : preview.error ? (
          <ErrorState error={preview.error} subject={path ?? undefined} />
        ) : data ? (
          <>
            {data.truncated && (
              <p className="rg-muted text-sm" role="status">
                Truncated: showing the first {formatBytes(data.limit)} of{' '}
                {formatBytes(data.size)}. Download the file to see all of it.
              </p>
            )}
            {data.content === '' ? (
              <p className="rg-muted text-sm">This file is empty.</p>
            ) : (
              <pre className="rg-code rg-file-preview" aria-label={`Contents of ${data.path}`}>
                {data.content}
              </pre>
            )}
          </>
        ) : null}

        <BaseModal.Footer>
          {path && (
            <DownloadLink className="rg-btn" href={localFileDownloadUrl(client, assetDigest, path)}>
              {path.split('/').pop()}
            </DownloadLink>
          )}
          <button type="button" className="rg-btn" onClick={onClose}>
            Close
          </button>
        </BaseModal.Footer>
      </div>
    </BaseModal>
  );
}
