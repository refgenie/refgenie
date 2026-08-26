import { useQuery } from '@tanstack/react-query';
import { useLocalApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import {
  listRemoteAssets,
  listRemoteGenomes,
  listRemoteServers,
} from '../../services/resources/remote';

export function useRemoteServers(options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useQuery({
    queryKey: qk.remoteServers(),
    queryFn: ({ signal }) => listRemoteServers(client, { signal }),
    enabled: options?.enabled ?? true,
  });
}

export function useRemoteGenomes(serverUrl: string | undefined, options?: { enabled?: boolean }) {
  const client = useLocalApiClient();
  return useQuery({
    queryKey: qk.remoteGenomes(serverUrl),
    queryFn: ({ signal }) => listRemoteGenomes(client, { server_url: serverUrl }, { signal }),
    enabled: options?.enabled ?? true,
  });
}

export function useRemoteAssets(
  serverUrl: string | undefined,
  genomeDigest: string | undefined,
  options?: { enabled?: boolean },
) {
  const client = useLocalApiClient();
  return useQuery({
    queryKey: qk.remoteAssets(serverUrl, genomeDigest ?? ''),
    queryFn: ({ signal }) =>
      listRemoteAssets(
        client,
        { genome_digest: genomeDigest as string, server_url: serverUrl },
        { signal },
      ),
    enabled: !!genomeDigest && (options?.enabled ?? true),
  });
}
