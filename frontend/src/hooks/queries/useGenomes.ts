import { useMemo } from 'react';
import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { fetchAllGenomes, getGenome, listGenomes } from '../../services/resources/genomes';
import { buildGenomeIndex } from '../../utils/genomes';
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
 * The whole-backend genome list as one cache entry, so /genomes, /species and
 * the build form share a single request. Pass a module-level `build` so the
 * memo holds across renders.
 */
export function useGenomeIndexResource<S>(build: (genomes: readonly GenomeResponse[]) => S) {
  const client = useApiClient();
  const query = useResource(qk.genomeIndex(), ({ signal }) => fetchAllGenomes(client, { signal }), {
    staleTime: INDEX_STALE_MS,
  });
  const data = useMemo(() => query.data && build(query.data), [query.data, build]);
  return { ...query, data };
}

/**
 * Digest -> genome for the whole backend, so any row carrying only a
 * `genome_digest` can render an alias.
 */
export function useGenomeIndex() {
  return useGenomeIndexResource(buildGenomeIndex);
}
