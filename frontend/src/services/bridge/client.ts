/**
 * The two ApiClients pointed at a connected local refgenie.
 *
 * Credential omission is a property of the CLIENT, not of each call site: a
 * wrapping fetch means no bridge request can forget it. Everything else —
 * the error envelope, the timeout, and the `X-Refgenie-Action` header on
 * `mutate` — comes from the shared `ApiClient` unchanged.
 *
 * The local refgenie IS a refgenie1 server, so its paths are the ones in
 * `services/contracts.ts`. There is no second path vocabulary.
 *
 * `ApiClient.mutate` sends `content-type: application/json` and
 * `X-Refgenie-Action`, which triggers a CORS preflight. The server allows
 * exactly those two headers plus POST, and sets
 * `Access-Control-Allow-Private-Network` when the bridge is on — so the
 * preflight, including Chrome's LNA leg, already succeeds. Do not add any
 * other header.
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
