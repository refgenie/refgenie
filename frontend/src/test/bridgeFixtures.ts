import { coerceCapabilities } from '../services/capabilities';
import { CAPABILITY_KEYS } from '../types/ui';
import type { PingResponse } from '../services/bridge/contract';

export const BRIDGE_ORIGIN = 'http://localhost:8080';

/** Every capability switched on, write ones included. */
export const ALL_CAPABILITIES = coerceCapabilities(
  Object.fromEntries(CAPABILITY_KEYS.map((key) => [key, true])),
);

/** A local refgenie's ping; defaults to a full bridge that can pull. */
export function makePing(overrides: Partial<PingResponse> = {}): PingResponse {
  return {
    service: 'refgenie',
    bridge_version: 1,
    mode: 'local',
    refgenie_version: '1.0.0',
    instance_id: 'instance-1',
    instance_label: 'my laptop refgenie',
    bridge_mode: 'full',
    action_header: 'X-Refgenie-Action',
    capabilities: coerceCapabilities({ pull: true }),
    bridge: { actions_cross_origin: true },
    ...overrides,
  };
}
