/**
 * The current browse scope's route prefix: `''` on the page's own backend,
 * `/local` inside the bridge's `/local/*` branch.
 */

import { useContext } from 'react';
import { BridgeScopeContext } from '../services/bridge/scope';

export function useRouteBase(): '' | '/local' {
  return useContext(BridgeScopeContext).routeBase;
}
