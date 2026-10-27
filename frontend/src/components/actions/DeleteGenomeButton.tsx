/**
 * Delete a genome.
 *
 * `genome.remove` cascades UNCONDITIONALLY — the genome, all its aliases, all
 * its asset groups and all its assets go, and there is no `force` flag gating
 * any of it. The modal therefore states the real consequence with the live
 * asset count and requires the user to type the primary alias.
 *
 * An unaliased genome (one just `init`ed, before `alias set`) has nothing but a
 * digest to identify it. The phrase to type is then the digest's first eight
 * characters, not all thirty-two: a confirmation gesture that can only be
 * completed by copy-paste confirms nothing. The full digest is shown in the
 * modal body so the user can still verify what they are about to destroy.
 */

import { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { deleteGenome } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ConfirmModal } from '../common/ConfirmModal';
import { DigestChip } from '../common/DigestChip';
import { formatDigest } from '../../utils/format';

/** Characters of the digest to type when a genome carries no alias. */
const DIGEST_CONFIRM_LENGTH = 8;

export interface DeleteGenomeButtonProps {
  digest: string;
  /** Empty when the genome has no alias; the digest then stands in. */
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

  // What we CALL the genome, and what the user has to TYPE, are not the same
  // thing once there is no alias.
  const phrase = primaryAlias || formatDigest(digest);
  const confirmation = primaryAlias || digest.slice(0, DIGEST_CONFIRM_LENGTH);

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
        confirmationText={confirmation}
        body={
          <>
            <p className="mb-2">
              This deletes the genome, {aliasCount} alias{aliasCount === 1 ? '' : 'es'}, every
              asset group, and all {assetCount} asset{assetCount === 1 ? '' : 's'} with their
              files on disk.
            </p>
            <p className="mb-2">
              The cascade is unconditional. There is no way to keep the assets.
            </p>
            <p className="flex items-center gap-2 flex-wrap">
              <span className="rg-muted">Genome digest</span>
              <DigestChip digest={digest} length={32} />
            </p>
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
