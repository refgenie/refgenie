/**
 * Submitting a pull.
 *
 * On the 202 this registers the job in the store immediately, using the id the
 * server just returned. That is not optimistic UI — the record is real; only
 * the label is provisional, and the first `status` frame replaces the whole
 * record. It exists so the console shows the job before the first SSE frame,
 * which on a cold connection can be a second away.
 *
 * There is no duplicate-submission error to handle: the manager coalesces and
 * answers 202 with `duplicate: true` and the SAME job id, so the hook focuses
 * the existing card instead of creating a second one.
 */

import { useCallback, useState } from 'react';
import { useLocalApiClient } from './useApiClient';
import { useToast } from './useToast';
import { pullAsset } from '../services/actions';
import { provisionalJob, useJobStore } from '../stores/jobStore';
import { ApiError } from '../services/http';
import type { JobRef } from '../services/contracts';

export interface PullInput {
  /**
   * The request model requires EXACTLY ONE genome reference: an alias in
   * `genome` or a digest in `genome_digest`. Sending both is a 422, and so is
   * sending neither.
   */
  genomeName?: string | null;
  genomeDigest?: string | null;
  assetGroupName: string;
  assetName?: string | null;
  serverUrl?: string | null;
  /**
   * ALWAYS explicit. The request model defaults it to `False` precisely so the
   * puller cannot reach its stdin prompt and hang a worker; a UI that relied on
   * that default would be one refactor away from the same hang.
   */
  force: boolean;
}

export function pullLabel(input: PullInput): string {
  const genome = input.genomeDigest ?? input.genomeName ?? '';
  const asset = input.assetName ? `:${input.assetName}` : '';
  const from = input.serverUrl ? ` from ${input.serverUrl}` : '';
  return `pull ${genome}/${input.assetGroupName}${asset}${from}`;
}

export interface UsePullActionResult {
  pull: (input: PullInput) => Promise<JobRef | null>;
  isPending: boolean;
}

export function usePullAction(): UsePullActionResult {
  const client = useLocalApiClient();
  const toast = useToast();
  const [isPending, setPending] = useState(false);

  const pull = useCallback(
    async (input: PullInput) => {
      setPending(true);
      try {
        // Exactly one genome reference: the digest wins when we have it, and
        // the other key is omitted entirely rather than sent as null.
        const ref = await pullAsset(client, {
          asset_group: input.assetGroupName,
          ...(input.genomeDigest
            ? { genome_digest: input.genomeDigest }
            : { genome: input.genomeName ?? '' }),
          asset: input.assetName ?? null,
          server_url: input.serverUrl ?? null,
          force: input.force,
        });

        const store = useJobStore.getState();
        if (ref.duplicate) {
          store.focusJob(ref.job_id);
          toast.info('That pull is already running.');
          return ref;
        }

        store.registerQueued(
          provisionalJob({
            id: ref.job_id,
            kind: ref.kind ?? 'pull',
            label: pullLabel(input),
            status: ref.status ?? 'queued',
            created_at: ref.created_at,
            target: {
              genome_digest: input.genomeDigest ?? null,
              genome_name: input.genomeName ?? null,
              asset_group_name: input.assetGroupName,
              asset_name: input.assetName ?? null,
            },
          }),
        );
        toast.success('Pull queued. Progress is in the job console.');
        return ref;
      } catch (error) {
        toast.error(
          error instanceof ApiError ? error.detail : 'Could not submit the pull.',
        );
        return null;
      } finally {
        setPending(false);
      }
    },
    [client, toast],
  );

  return { pull, isPending };
}
