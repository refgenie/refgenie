/**
 * Delete an asset. Synchronous, not a job: `asset.remove_by_digest` is a DB
 * delete plus a directory unlink, which finishes inside the request.
 */

import { useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { deleteAsset } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ConfirmModal } from '../common/ConfirmModal';
import type { SeekKeyResponse } from '../../types/api';

export interface DeleteAssetButtonProps {
  digest: string;
  registryPath: string;
  seekKeys?: SeekKeyResponse[] | null;
  onDeleted?: () => void;
  compact?: boolean;
}

export function DeleteAssetButton({
  digest,
  registryPath,
  seekKeys,
  onDeleted,
  compact = false,
}: DeleteAssetButtonProps) {
  const canDelete = useCapability('delete');
  const client = useLocalApiClient();
  const toast = useToast();
  const [open, setOpen] = useState(false);

  if (!canDelete) return null;

  return (
    <>
      <button
        type="button"
        className={compact ? 'rg-btn rg-btn--sm rg-btn--danger' : 'rg-btn rg-btn--danger'}
        onClick={() => setOpen(true)}
      >
        Delete
      </button>

      <ConfirmModal
        isOpen={open}
        onClose={() => setOpen(false)}
        title="Delete asset"
        destructive
        confirmLabel="Delete asset"
        body={
          <>
            <p className="mb-2">
              Delete <code className="rg-code rg-code--inline">{registryPath}</code>? Its files
              are removed from disk.
            </p>
            {seekKeys && seekKeys.length > 0 && (
              <>
                <p className="mb-1">These seek keys will disappear:</p>
                <ul className="rg-list">
                  {seekKeys.map((key) => (
                    <li key={key.name}>
                      <code className="rg-code rg-code--inline">{key.name}</code>
                    </li>
                  ))}
                </ul>
              </>
            )}
          </>
        }
        onConfirm={async () => {
          await deleteAsset(client, digest);
          invalidate(INVALIDATION_FOR_ACTION['asset.delete']);
          toast.success(`Deleted ${registryPath}.`);
          onDeleted?.();
        }}
      />
    </>
  );
}
