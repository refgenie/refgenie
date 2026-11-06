/**
 * ApiClients for a connected local refgenie, with credentials always omitted.
 * The server's CORS policy allows only `content-type` and `X-Refgenie-Action`,
 * so do not add other headers.
 */

import { ApiClient } from '../http';

const omitCredentials: typeof fetch = (input, init) =>
  globalThis.fetch(input, { ...init, credentials: 'omit' });

export interface BridgeClients {
  /** Shared read API on the local instance. */
  read: ApiClient;
  /** Local-only surface: `POST /v1/actions/pull` and `GET /v1/jobs/{id}`. */
  local: ApiClient;
}

export function createBridgeClients(baseUrl: string): BridgeClients {
  return {
    read: new ApiClient({ baseUrl: `${baseUrl}/v4`, fetchImpl: omitCredentials }),
    local: new ApiClient({ baseUrl: `${baseUrl}/v1`, fetchImpl: omitCredentials }),
  };
}
