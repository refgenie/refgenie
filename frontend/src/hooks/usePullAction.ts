/** Submits a pull and puts its job on the console. */

import { useCallback } from 'react';
import { useLocalApiClient } from './useApiClient';
import { useToast } from './useToast';
import { pullAsset } from '../services/actions';
import { trackJob } from '../stores/jobStore';
import { useAction } from './useAction';

export interface PullInput {
  /** Give a name or a digest; the request takes exactly one genome reference. */
  genomeName?: string | null;
  genomeDigest?: string | null;
  assetGroupName: string;
  assetName?: string | null;
  serverUrl?: string | null;
  /** Always explicit, so the puller can never hang on a stdin prompt. */
  force: boolean;
}

function pullLabel(input: PullInput): string {
  const genome = input.genomeDigest ?? input.genomeName ?? '';
  const asset = input.assetName ? `:${input.assetName}` : '';
  const from = input.serverUrl ? ` from ${input.serverUrl}` : '';
  return `pull ${genome}/${input.assetGroupName}${asset}${from}`;
}

export interface UsePullActionResult {
  /** Resolves true when the pull was submitted. */
  pull: (input: PullInput) => Promise<boolean>;
  isPending: boolean;
}

export function usePullAction(): UsePullActionResult {
  const client = useLocalApiClient();
  const toast = useToast();
  const { busy: isPending, run } = useAction(toast.error, 'Could not submit the pull.');

  const pull = useCallback(
    (input: PullInput) =>
      run(async () => {
        // Exactly one genome reference: prefer the digest, and omit the other key.
        const ref = await pullAsset(client, {
          asset_group: input.assetGroupName,
          ...(input.genomeDigest
            ? { genome_digest: input.genomeDigest }
            : { genome: input.genomeName ?? '' }),
          asset: input.assetName ?? null,
          server_url: input.serverUrl ?? null,
          force: input.force,
        });
        trackJob(
          ref,
          {
            noun: 'pull',
            label: pullLabel(input),
            target: {
              genome_digest: input.genomeDigest ?? null,
              genome_name: input.genomeName ?? null,
              asset_group_name: input.assetGroupName,
              asset_name: input.assetName ?? null,
            },
          },
          toast,
        );
      }),
    [client, toast, run],
  );

  return { pull, isPending };
}
