import { createContext } from 'react';

/**
 * Which backend the surrounding subtree is browsing, expressed as the route
 * prefix its internal links must keep. `''` is the page's own backend;
 * `/local` is a connected local refgenie reached over the bridge.
 */
export interface BridgeScope {
  routeBase: '' | '/local';
}

export const BridgeScopeContext = createContext<BridgeScope>({ routeBase: '' });
