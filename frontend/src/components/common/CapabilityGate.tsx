import type { ReactNode } from 'react';
import { useCapability } from '../../hooks/useCapability';
import type { CapabilityKey } from '../../types/ui';

export interface CapabilityGateProps {
  cap: CapabilityKey;
  children: ReactNode;
}

/**
 * Renders children only when the capability is on, and `null` otherwise. There
 * are no disabled-but-visible buttons in v1.
 */
export function CapabilityGate({ cap, children }: CapabilityGateProps) {
  return useCapability(cap) ? <>{children}</> : null;
}
