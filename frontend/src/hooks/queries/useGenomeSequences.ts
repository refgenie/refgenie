import { useResource } from '../useResource';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { NO_RETRY } from '../../services/resourceCache';
import { compareGenomes, getGenomeSequences } from '../../services/resources/sequences';

/** A 404 means "no sequence data here", a normal answer, so never retry. */
export function useGenomeSequences(digest: string | undefined) {
  const client = useApiClient();
  return useResource(
    qk.genomeSequences(digest ?? ''),
    ({ signal }) => getGenomeSequences(client, digest as string, { signal }),
    { enabled: !!digest, retry: NO_RETRY },
  );
}

/** Comparing a genome with itself is allowed. */
export function useGenomeComparison(a: string | undefined, b: string | undefined) {
  const client = useApiClient();
  return useResource(
    qk.genomeComparison(a ?? '', b ?? ''),
    ({ signal }) => compareGenomes(client, a as string, b as string, { signal }),
    { enabled: !!a && !!b, retry: NO_RETRY },
  );
}
