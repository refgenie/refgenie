/**
 * Remote browse. Local mode only, gated on `capabilities.remote_browse`.
 *
 * These live on the local-only `/v1` surface, so they take the *local* client
 * (`useLocalApiClient()`), not the shared read client.
 */

import type { ApiClient, RequestInitLite } from '../http';
import type {
  RemoteAsset,
  RemoteGenome,
  RemoteServersResponse,
} from '../../types/api';

export const listRemoteServers = (c: ApiClient, init?: RequestInitLite) =>
  c.get<RemoteServersResponse>('/remote/servers', undefined, init);

export const listRemoteGenomes = (
  c: ApiClient,
  p: { server_url?: string } = {},
  init?: RequestInitLite,
) => c.get<RemoteGenome[]>('/remote/genomes', { server_url: p.server_url }, init);

export const listRemoteAssets = (
  c: ApiClient,
  p: { genome_digest: string; server_url?: string },
  init?: RequestInitLite,
) =>
  c.get<RemoteAsset[]>(
    '/remote/assets',
    { genome_digest: p.genome_digest, server_url: p.server_url },
    init,
  );
