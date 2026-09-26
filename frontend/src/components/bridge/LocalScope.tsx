/**
 * The `/local/*` route branch's scope: the existing browse components render
 * against a connected local refgenie unchanged, because they are
 * client-parameterized through AppContext and never special-cased.
 *
 * Three things are swapped, and each one matters:
 *
 *  - AppContext, so `useApiClient()` and `useUiConfig()` resolve to the bridge
 *    target. The capability set is masked READ-ONLY: the ping's flags say what
 *    that instance can do, but the only state change this app ever makes
 *    cross-origin is `POST /v1/actions/pull`, which lives on the remote branch.
 *    Rendering Delete or Build here would offer an action the server is
 *    designed to refuse.
 *  - A separate resource cache. Cache keys in `services/queryKeys.ts` are NOT
 *    prefixed with a base URL, so sharing the root cache would collide
 *    `/local/genomes` with `/genomes` — the local scope would show remote data
 *    and the invalidation bus would cross-invalidate. The previous SPA solved
 *    this by prefixing every key with `api.baseUrl`; a second cache is the
 *    smaller change and isolates completely.
 *  - BridgeScopeContext, so internal links stay under `/local`.
 *
 * Local data is always visually segregated (banner + `/local` routes) and never
 * merged into remote result sets without a badge: a rogue process on the bridge
 * port could feed fabricated data, and presentation is the only available
 * mitigation (threat T7).
 */

import { useEffect, useMemo, useState } from 'react';
import type { ReactNode } from 'react';
import { createAppResourceCache } from '../../app/resourceCache';
import { ResourceCacheContext } from '../../services/resourceCacheContext';
import { AppContext } from '../../services/clients';
import { createBridgeClients } from '../../services/bridge/client';
import { BridgeScopeContext } from '../../services/bridge/scope';
import { noCapabilities } from '../../services/capabilities';
import { ConnectDialog } from './ConnectDialog';
import { useBridge } from '../../hooks/useBridge';
import type { UiConfig } from '../../types/ui';

const LOCAL_SCOPE = { routeBase: '/local' } as const;

export interface LocalScopeProps {
  children: ReactNode;
}

export function LocalScope({ children }: LocalScopeProps) {
  const bridge = useBridge();
  const [dialogOpen, setDialogOpen] = useState(false);
  const baseUrl = bridge.baseUrl;
  const ping = bridge.ping;

  const clients = useMemo(
    () => (baseUrl ? createBridgeClients(baseUrl) : null),
    [baseUrl],
  );

  // A second cache, torn down when the connection target changes: cache keys
  // are not base-URL-prefixed, so this is what keeps `/local` isolated.
  const resourceCache = useMemo(() => {
    // Keyed on the target: reconnecting to a different port must not reuse a
    // cache populated from the previous one.
    void baseUrl;
    return createAppResourceCache();
  }, [baseUrl]);
  useEffect(() => () => resourceCache.clear(), [resourceCache]);

  const localConfig: UiConfig = useMemo(
    () => ({
      mode: 'local',
      api_base: `${baseUrl ?? ''}/v4`,
      root_path: '',
      refgenie_version: ping?.refgenie_version ?? 'unknown',
      service_name: ping?.instance_label ?? 'local refgenie',
      capabilities: {
        // Every command flag stays false. That is the mask, and it is
        // intentional: do not "fix" it later by spreading the ping wholesale.
        ...noCapabilities(),
        downloads: ping?.capabilities.downloads === true,
        archives: ping?.capabilities.archives === true,
        seqcol: ping?.capabilities.seqcol === true,
        drs: ping?.capabilities.drs === true,
      },
      web_ui: null,
      degraded: false,
    }),
    [baseUrl, ping],
  );

  const appValue = useMemo(
    () =>
      clients ? { config: localConfig, client: clients.read, localClient: clients.local } : null,
    [clients, localConfig],
  );

  if (bridge.status !== 'connected' || !appValue) {
    return (
      <div className="rg-card">
        <div className="rg-card__body flex flex-col gap-3">
          <h2 className="text-xl font-semibold">Not connected to a local refgenie</h2>
          <p className="rg-muted text-sm">
            {bridge.status === 'connecting'
              ? 'Connecting…'
              : 'Connect to your running refgenie dash to browse your local genomes here.'}
          </p>
          {bridge.status !== 'connecting' && (
            <span>
              <button
                type="button"
                className="rg-btn rg-btn--primary"
                onClick={() => setDialogOpen(true)}
              >
                Connect
              </button>
            </span>
          )}
        </div>
        <ConnectDialog open={dialogOpen} onClose={() => setDialogOpen(false)} />
      </div>
    );
  }

  return (
    <ResourceCacheContext.Provider value={resourceCache}>
      <AppContext.Provider value={appValue}>
        <BridgeScopeContext.Provider value={LOCAL_SCOPE}>
          <div className="flex flex-col gap-6">
            <div className="rg-banner rg-banner--info" role="status">
              <span>
                Browsing your <strong>local refgenie</strong> ({ping?.instance_label}) —
                this data lives on your computer.
              </span>
            </div>
            {children}
          </div>
        </BridgeScopeContext.Provider>
      </AppContext.Provider>
    </ResourceCacheContext.Provider>
  );
}
