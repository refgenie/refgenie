/**
 * The `/ping` handshake of the localhost bridge.
 *
 * `bridge_version` changes only on breaking shape changes. Gate features on
 * `capabilities` (a missing flag is false), never on version numbers.
 */

import { coerceCapabilities } from '../capabilities';
import type { Capabilities } from '../../types/ui';

export const SUPPORTED_BRIDGE_VERSIONS: ReadonlySet<number> = new Set([1]);

export type PingResponse = {
  service: 'refgenie';
  bridge_version: number;
  mode: string;
  refgenie_version: string;
  api_version?: string;
  instance_id: string;
  instance_label: string;
  bridge_mode: 'off' | 'read' | 'full' | string;
  action_header: string;
  capabilities: Capabilities;
  bridge: { actions_cross_origin?: boolean };
  /** Absent when the server's count queries fail. */
  counts?: { genomes: number; assets: number };
};

export type PingValidation =
  | { ok: true; ping: PingResponse }
  | { ok: false; reason: 'wrong-service' | 'unsupported-version' | 'malformed' };

const isRecord = (value: unknown): value is Record<string, unknown> =>
  typeof value === 'object' && value !== null && !Array.isArray(value);

/**
 * Validate an untrusted `/ping` payload; any program could be on the port. This
 * guards against port collisions, not attacks: a page cannot authenticate a
 * loopback peer, so local data is always shown badged and apart.
 */
export const validatePing = (data: unknown): PingValidation => {
  if (!isRecord(data)) return { ok: false, reason: 'malformed' };
  if (data.service !== 'refgenie') return { ok: false, reason: 'wrong-service' };
  if (typeof data.bridge_version !== 'number') {
    return { ok: false, reason: 'malformed' };
  }
  if (!SUPPORTED_BRIDGE_VERSIONS.has(data.bridge_version)) {
    return { ok: false, reason: 'unsupported-version' };
  }
  for (const field of [
    'mode',
    'instance_id',
    'instance_label',
    'bridge_mode',
    'action_header',
  ]) {
    if (typeof data[field] !== 'string') return { ok: false, reason: 'malformed' };
  }
  if (!isRecord(data.capabilities) || !isRecord(data.bridge)) {
    return { ok: false, reason: 'malformed' };
  }
  const capabilities = coerceCapabilities(
    data.capabilities as Partial<Record<string, boolean>>,
  );
  return {
    ok: true,
    ping: {
      service: 'refgenie',
      bridge_version: data.bridge_version,
      mode: data.mode as string,
      refgenie_version:
        typeof data.refgenie_version === 'string' ? data.refgenie_version : 'unknown',
      api_version: typeof data.api_version === 'string' ? data.api_version : undefined,
      instance_id: data.instance_id as string,
      instance_label: data.instance_label as string,
      bridge_mode: data.bridge_mode as string,
      action_header: data.action_header as string,
      capabilities,
      bridge: {
        actions_cross_origin:
          (data.bridge as Record<string, unknown>).actions_cross_origin === true,
      },
      counts:
        isRecord(data.counts) &&
        typeof data.counts.genomes === 'number' &&
        typeof data.counts.assets === 'number'
          ? { genomes: data.counts.genomes, assets: data.counts.assets }
          : undefined,
    },
  };
};
