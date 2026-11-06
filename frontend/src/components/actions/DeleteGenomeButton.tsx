/**
 * Delete a genome. This always removes its aliases, asset groups and assets
 * too, so the modal shows the asset count and asks the user to type the alias
 * (or the first 8 digest characters if the genome has no alias).
 */

import { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { deleteGenome } from '../../services/actions';
import { useInvalidate } from '../../hooks/useResource';
import { staleAfter } from '../../services/queryKeys';
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
  const invalidate = useInvalidate();
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
          invalidate(staleAfter['genome.delete']);
          toast.success(`Deleted ${phrase}.`);
          navigate('/genomes');
        }}
      />
    </>
  );
}
