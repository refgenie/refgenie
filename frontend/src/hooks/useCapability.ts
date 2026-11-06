import { useUiConfig } from './useUiConfig';
import type { CapabilityKey } from '../types/ui';

/** Whether a capability is on. Gate commands on capabilities, never on the mode string. */
export function useCapability(cap: CapabilityKey): boolean {
  return useUiConfig().capabilities[cap] === true;
}
