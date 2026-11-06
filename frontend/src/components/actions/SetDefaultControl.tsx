/**
 * Picks the group's default asset. The radio changes only after the server
 * confirms, not optimistically.
 */

import { useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { setDefaultAsset } from '../../services/actions';
import { useInvalidate } from '../../hooks/useResource';
import { staleAfter } from '../../services/queryKeys';
import { useAction } from '../../hooks/useAction';

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
  const invalidate = useInvalidate();
  const { busy, run } = useAction(toast.error, 'Could not set the default asset.');
  const [confirmed, setConfirmed] = useState(isDefault);

  if (!canCurate) return null;

  const id = `default-${assetGroupName}-${assetName}`;

  const submit = () =>
    run(async () => {
      await setDefaultAsset(client, {
        genome_digest: genomeDigest,
        asset_group: assetGroupName,
        asset: assetName,
      });
      setConfirmed(true);
      invalidate(staleAfter['asset.set_default']);
      toast.success(`${assetGroupName}:${assetName} is now the default.`);
    });

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
