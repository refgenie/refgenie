import { createContext } from 'react';

/**
 * Route prefix for the backend a subtree browses: `''` for this page's backend,
 * `/local` for a connected local refgenie.
 */
export interface BridgeScope {
  routeBase: '' | '/local';
}

export const BridgeScopeContext = createContext<BridgeScope>({ routeBase: '' });
