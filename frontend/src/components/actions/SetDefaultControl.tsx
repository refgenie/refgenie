/**
 * Which asset of a group is the default one a bare registry path resolves to.
 *
 * The radio flips only AFTER a 200: this touches the database, and an
 * optimistic flip that the server then rejected would leave the list lying
 * about which asset is default.
 */

import { useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { setDefaultAsset } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ApiError } from '../../services/http';

export interface SetDefaultControlProps {
  genomeDigest: string;
  assetGroupName: string;
  assetName: string;
  isDefault: boolean;
}

export function SetDefaultControl({
  genomeDigest,
  assetGroupName,
  assetName,
  isDefault,
}: SetDefaultControlProps) {
  // Asset curation travels with the build capability in the shared key set.
  const canCurate = useCapability('build');
  const client = useLocalApiClient();
  const toast = useToast();
  const [busy, setBusy] = useState(false);
  const [confirmed, setConfirmed] = useState(isDefault);

  if (!canCurate) return null;

  const id = `default-${assetGroupName}-${assetName}`;

  const submit = async () => {
    setBusy(true);
    try {
      await setDefaultAsset(client, {
        genome_digest: genomeDigest,
        asset_group: assetGroupName,
        asset: assetName,
      });
      setConfirmed(true);
      invalidate(INVALIDATION_FOR_ACTION['asset.set_default']);
      toast.success(`${assetGroupName}:${assetName} is now the default.`);
    } catch (error) {
      toast.error(
        error instanceof ApiError ? error.detail : 'Could not set the default asset.',
      );
    } finally {
      setBusy(false);
    }
  };

  return (
    <label className="rg-check rg-check--inline" htmlFor={id}>
      <input
        id={id}
        type="radio"
        name={`default-asset-${assetGroupName}`}
        checked={confirmed}
        disabled={busy}
        onChange={() => void submit()}
      />
      <span className="sr-only">Make {assetName} the default asset</span>
    </label>
  );
}
