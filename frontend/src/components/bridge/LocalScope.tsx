/**
 * Scope for the `/local/*` routes: the normal browse pages render against a
 * connected local refgenie by swapping AppContext (read-only capabilities),
 * the resource cache (keys are not prefixed by base URL, so they would collide)
 * and BridgeScopeContext (links stay under `/local`). A banner marks the data as
 * local, since a rogue process on the bridge port could serve fake data.
 */

import { useEffect, useMemo, useState } from 'react';
import type { ReactNode } from 'react';
import { NavLink, useLocation } from 'react-router-dom';
import { createAppResourceCache } from '../../app/resourceCache';
import { ResourceCacheContext } from '../../hooks/useResource';
import { AppContext } from '../../services/clients';
import { createBridgeClients } from '../../services/bridge/client';
import { BridgeScopeContext } from '../../services/bridge/scope';
import { noCapabilities } from '../../services/capabilities';
import { ConnectDialog } from './ConnectDialog';
import { useBridge } from '../../hooks/useBridge';
import { cn } from '../../utils/cn';
import { NAV } from '../../app/navigation';
import type { UiConfig } from '../../types/ui';

const LOCAL_SCOPE = { routeBase: '/local' } as const;

/** Local navigation; the top bar links point at the public catalog. */
const LOCAL_NAV = NAV.filter((entry) => entry.local).map((entry) => ({
  to: `/local${entry.to}`,
  label: entry.label,
}));

/** Bare `/local` renders the genomes list, so its nav entry is the current one. */
function isCurrent(to: string, pathname: string, isActive: boolean): boolean {
  return isActive || (to === '/local/genomes' && /^\/local\/?$/.test(pathname));
}

export interface LocalScopeProps {
  children: ReactNode;
}

export function LocalScope({ children }: LocalScopeProps) {
  const bridge = useBridge();
  const { pathname } = useLocation();
  const [dialogOpen, setDialogOpen] = useState(false);
  const baseUrl = bridge.baseUrl;
  const ping = bridge.ping;

  const clients = useMemo(
    () => (baseUrl ? createBridgeClients(baseUrl) : null),
    [baseUrl],
  );

  // A new cache per target, so reconnecting to another port starts empty.
  const resourceCache = useMemo(() => {
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
        // Command flags stay false: the only cross-origin write is pull, on the remote branch.
        ...noCapabilities(),
        // `jobs` gates only the read-only job history page here.
        jobs: ping?.capabilities.jobs === true,
        downloads: ping?.capabilities.downloads === true,
        archives: ping?.capabilities.archives === true,
        seqcol: ping?.capabilities.seqcol === true,
        drs: ping?.capabilities.drs === true,
      },
      web_ui: null,
      collections: [],
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
      <div>
        <section className="flex flex-col gap-3">
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
        </section>
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
                Browsing your <strong>local refgenie</strong> ({ping?.instance_label}). This data lives on your computer.
              </span>
            </div>
            <nav aria-label="Local refgenie">
              <ul className="rg-nav__list">
                {LOCAL_NAV.map((entry) => (
                  <li key={entry.to}>
                    <NavLink
                      to={entry.to}
                      className={({ isActive }) =>
                        cn(
                          'rg-nav__item',
                          isCurrent(entry.to, pathname, isActive) && 'rg-nav__item--active',
                        )
                      }
                    >
                      {entry.label}
                    </NavLink>
                  </li>
                ))}
              </ul>
            </nav>
            {children}
          </div>
        </BridgeScopeContext.Provider>
      </AppContext.Provider>
    </ResourceCacheContext.Provider>
  );
}
