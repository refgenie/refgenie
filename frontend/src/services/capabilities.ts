/**
 * Capability parsing for `/service-info` and `/ping`. Only `true` counts as on;
 * keys outside `CAPABILITY_KEYS` are dropped.
 */

import { CAPABILITY_KEYS } from '../types/ui';
import type { Capabilities } from '../types/ui';

/** Every capability false. */
export function noCapabilities(): Capabilities {
  return Object.fromEntries(
    CAPABILITY_KEYS.map((k) => [k, false]),
  ) as unknown as Capabilities;
}

export function coerceCapabilities(
  raw: Partial<Record<string, boolean>> | undefined,
): Capabilities {
  const out = noCapabilities();
  if (!raw) return out;
  for (const key of CAPABILITY_KEYS) {
    out[key] = raw[key] === true;
  }
  return out;
}
