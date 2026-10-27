/**
 * The React context that carries the active `ResourceCache`.
 *
 * Kept beside `clients.ts` and shaped the same way: a bare context, no provider
 * component, so the three places that scope a cache use
 * `<ResourceCacheContext.Provider>` directly and the hook file stays free of
 * components.
 *
 * Null outside a provider, which is a programming error — see
 * `hooks/useResource.ts`.
 */

import { createContext } from 'react';
import type { ResourceCache } from './resourceCache';

export const ResourceCacheContext = createContext<ResourceCache | null>(null);
