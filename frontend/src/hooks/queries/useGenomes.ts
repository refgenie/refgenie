import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { fetchAllGenomes, getGenome, listGenomes } from '../../services/resources/genomes';
import { buildGenomeIndex } from '../../utils/genomes';
import type { UseResourceResult } from '../useResource';
import type { ListGenomesParams } from '../../services/resources/genomes';
import type { GenomeResponse } from '../../types/api';

export function useGenomes(p: ListGenomesParams, options?: { enabled?: boolean }) {
  const client = useApiClient();
  return useResource(qk.genomes(p), ({ signal }) => listGenomes(client, p, { signal }), {
    enabled: options?.enabled ?? true,
  });
}

export function useGenome(digest: string | undefined) {
  const client = useApiClient();
  return useResource(
    qk.genome(digest ?? ''),
    ({ signal }) => getGenome(client, digest as string, { signal }),
    { enabled: !!digest },
  );
}

const INDEX_STALE_MS = 5 * 60 * 1000;

/**
 * The whole-backend genome list, as ONE cache entry.
 *
 * `useGenomeIndex` and `useSpeciesIndex` are two `select`s over this. The cache
 * stores the raw array under `qk.genomeIndex()` and memoizes each selector
 * separately (keyed on the selector's own identity, which is why both must be
 * module-level functions), so opening /species after /genomes costs no request
 * and builds neither Map twice. One request per session; see `fetchAllGenomes`.
 *
 * Both public hooks below go through here rather than calling `useResource`
 * themselves, so there is exactly one place the key, the fetcher and the stale
 * time are written down.
 */
export function useGenomeIndexResource<S>(
  select: (genomes: readonly GenomeResponse[]) => S,
): UseResourceResult<S> {
  const client = useApiClient();
  return useResource(qk.genomeIndex(), ({ signal }) => fetchAllGenomes(client, { signal }), {
    select,
    staleTime: INDEX_STALE_MS,
  });
}

/**
 * Digest -> genome for the whole backend, so any row carrying only a
 * `genome_digest` can render an alias.
 */
export function useGenomeIndex() {
  return useGenomeIndexResource(buildGenomeIndex);
}
