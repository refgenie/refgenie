import { useResource } from '../useResource';
import { useLocalApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import {
  listDataChannels,
  listRemoteAssets,
  listRemoteGenomes,
  listRemoteServers,
} from '../../services/resources/remote';

export function useRemoteServers(options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useResource(qk.remoteServers(), ({ signal }) => listRemoteServers(client, { signal }), {
    enabled: options?.enabled ?? true,
  });
}

export function useDataChannels(options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useResource(qk.dataChannels(), ({ signal }) => listDataChannels(client, { signal }), {
    enabled: options?.enabled ?? true,
  });
}

export function useRemoteGenomes(serverUrl: string | undefined, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useResource(
    qk.remoteGenomes(serverUrl),
    ({ signal }) => listRemoteGenomes(client, { server_url: serverUrl }, { signal }),
    { enabled: options?.enabled ?? true },
  );
}

export function useRemoteAssets(
  serverUrl: string | undefined,
  genomeDigest: string | undefined,
  options?: { enabled?: boolean },
) {
  const client = useLocalApiClient();
  return useResource(
    qk.remoteAssets(serverUrl, genomeDigest ?? ''),
    ({ signal }) =>
      listRemoteAssets(
        client,
        { genome_digest: genomeDigest as string, server_url: serverUrl },
        { signal },
      ),
    { enabled: !!genomeDigest && (options?.enabled ?? true) },
  );
}
