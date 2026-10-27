/**
 * Capability coercion, shared by the two documents that carry capability
 * flags: `/service-info` (the page's own backend) and `/ping` (a connected
 * local refgenie). Keeping one implementation is what stops the bridge's
 * copy from drifting away from the bootstrap contract.
 *
 * A missing key reads false. A non-boolean value reads false — never truthy.
 *
 * Keys outside `CAPABILITY_KEYS` are dropped rather than carried: an
 * undeclared flag has no consumer, and a capability is gateable only once it
 * is declared in the shared key set.
 */

import { CAPABILITY_KEYS } from '../types/ui';
import type { Capabilities } from '../types/ui';

/** Every capability false. A missing key always reads as false. */
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
