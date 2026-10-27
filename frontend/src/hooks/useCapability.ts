import { useUiConfig } from './useUiConfig';
import type { CapabilityKey } from '../types/ui';

/**
 * Every command affordance is gated on a capability flag, never on the mode
 * string. That is what lets a future authenticated server flip individual
 * flags without a frontend change.
 */
export function useCapability(cap: CapabilityKey): boolean {
  return useUiConfig().capabilities[cap] === true;
}
