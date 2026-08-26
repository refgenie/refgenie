/**
 * Delete a genome.
 *
 * `genome.remove` cascades UNCONDITIONALLY — the genome, all its aliases, all
 * its asset groups and all its assets go, and there is no `force` flag gating
 * any of it. The modal therefore states the real consequence with the live
 * asset count and requires the user to type the primary alias.
 */

import { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { deleteGenome } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ConfirmModal } from '../common/ConfirmModal';

export interface DeleteGenomeButtonProps {
  digest: string;
  /** Typed by the user to enable the confirm button. Falls back to the digest. */
  primaryAlias: string;
  aliasCount: number;
  assetCount: number;
}

export function DeleteGenomeButton({
  digest,
  primaryAlias,
  aliasCount,
  assetCount,
}: DeleteGenomeButtonProps) {
  const canDelete = useCapability('delete');
  const client = useLocalApiClient();
  const toast = useToast();
  const navigate = useNavigate();
  const [open, setOpen] = useState(false);

  if (!canDelete) return null;

  const phrase = primaryAlias || digest;

  return (
    <>
      <button type="button" className="rg-btn rg-btn--danger" onClick={() => setOpen(true)}>
        Delete genome
      </button>

      <ConfirmModal
        isOpen={open}
        onClose={() => setOpen(false)}
        title="Delete genome"
        destructive
        confirmLabel="Delete genome and all its assets"
        confirmationText={phrase}
        body={
          <>
            <p className="mb-2">
              This deletes the genome, {aliasCount} alias{aliasCount === 1 ? '' : 'es'}, every
              asset group, and all {assetCount} asset{assetCount === 1 ? '' : 's'} with their
              files on disk.
            </p>
            <p>The cascade is unconditional. There is no way to keep the assets.</p>
          </>
        }
        onConfirm={async () => {
          await deleteGenome(client, digest);
          invalidate(INVALIDATION_FOR_ACTION['genome.delete']);
          toast.success(`Deleted ${phrase}.`);
          navigate('/genomes');
        }}
      />
    </>
  );
}
