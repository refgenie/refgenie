import { validatePing } from './contract';
import type { PingResponse } from './contract';

/**
 * Probe one port, never a scan: scanning a visitor's ports looks like an attack,
 * and browsers ask permission per target. This is the only `fetch` outside
 * `ApiClient`, because it must stay a CORS simple request (no custom headers,
 * no credentials).
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
   * The fetch failed. Browsers report mixed-content blocks, permission denials,
   * CORS errors and refused connections the same way, so the UI lists all causes.
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
    // No credentials or custom headers, so no CORS preflight is needed.
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
