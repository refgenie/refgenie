import { validatePing } from './contract';
import type { PingResponse } from './contract';

/**
 * One probe of exactly one port — never a scan (D2): a public page
 * enumerating ports on the visitor's machine IS the attack this feature must
 * be distinguishable from, and Chrome's Local Network Access model multiplies
 * permission friction per distinct target.
 *
 * This is the ONE place in the app that calls `fetch` outside `ApiClient`, and
 * deliberately so. The probe's correctness is its request shape: a CORS
 * *simple* request with no custom headers and `credentials: 'omit'`, on a
 * 1.5 s budget. Routing it through the shared client would make that property
 * depend on internals that exist to serve a different purpose. Every other
 * bridge request DOES go through `ApiClient` — see `client.ts`.
 * `eslint.config.js` exempts `src/services/**` from the no-fetch rule.
 */

export const PROBE_TIMEOUT_MS = 1500;

export type ProbeResult =
  | { kind: 'connected'; ping: PingResponse; baseUrl: string }
  /** Something answered but it is not a refgenie we can talk to. */
  | {
      kind: 'unsupported';
      reason: 'wrong-service' | 'unsupported-version';
    }
  /**
   * The fetch rejected. Browsers collapse mixed-content blocks, LNA
   * permission denials, CORS rejections and connection-refused into an opaque
   * `TypeError: Failed to fetch` — the UI must present ALL real causes, not
   * pretend to know which one occurred.
   */
  | { kind: 'unreachable' };

export const localBaseUrl = (port: number): string => `http://localhost:${port}`;

export const probeLocal = async (
  port: number,
  fetchImpl: typeof fetch = fetch,
): Promise<ProbeResult> => {
  const baseUrl = localBaseUrl(port);
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), PROBE_TIMEOUT_MS);
  try {
    // credentials omitted and NO custom headers: the probe must stay a CORS
    // *simple* request so detection never depends on a preflight succeeding.
    const response = await fetchImpl(`${baseUrl}/ping`, {
      credentials: 'omit',
      signal: controller.signal,
    });
    if (!response.ok) return { kind: 'unreachable' };
    const validation = validatePing(await response.json());
    if (!validation.ok) {
      return {
        kind: 'unsupported',
        reason: validation.reason === 'malformed' ? 'wrong-service' : validation.reason,
      };
    }
    return { kind: 'connected', ping: validation.ping, baseUrl };
  } catch {
    return { kind: 'unreachable' };
  } finally {
    clearTimeout(timer);
  }
};
