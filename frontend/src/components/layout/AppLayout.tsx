import { useState } from 'react';
import { Outlet } from 'react-router-dom';
import { NavBar } from './NavBar';
import { Footer } from './Footer';
import { BridgeControl } from '../bridge/BridgeControl';
import { JobConsole } from '../jobs/JobConsole';
import { cn } from '../../utils/cn';
import { useUiConfig } from '../../hooks/useUiConfig';
import { useCapability } from '../../hooks/useCapability';
import { useInvalidationBridge } from '../../hooks/useInvalidationBridge';
import { useJobEvents } from '../../hooks/useJobEvents';
import { useBridgeAutoConnect } from '../../hooks/useBridgeAutoConnect';

export function AppLayout() {
  const config = useUiConfig();
  const hasJobs = useCapability('jobs');
  const [bannerDismissed, setBannerDismissed] = useState(false);

  // The single event transport and the single query-invalidation listener are
  // both mounted here, once. Nothing else in the tree opens a connection.
  useInvalidationBridge();
  useJobEvents({ enabled: hasJobs });
  // The bridge's one remembered-connection re-probe. A local dash IS the local
  // refgenie, so nothing bridge-related runs there at all.
  useBridgeAutoConnect(config.mode !== 'local');

  return (
    <div className={cn('rg-layout', hasJobs && 'rg-layout--docked')}>
      <header className="rg-layout__header">
        <div className="rg-layout__bar">
          <NavBar />
          {config.mode !== 'local' && <BridgeControl />}
        </div>
      </header>
      <div className="rg-layout__body">
        <main className="rg-layout__main">
          {/* Outside the content column, so a page's full-bleed head is always
              the first thing under the nav and can sit flush against it. */}
          {config.degraded && !bannerDismissed && (
            <div className="rg-banner rg-banner--warning rg-banner--strip" role="status">
              <span>
                Could not read <code className="rg-code rg-code--inline">/service-info</code>, so
                this instance&rsquo;s capabilities are unknown. The interface is read-only.
              </span>
              <button
                type="button"
                className="rg-btn rg-btn--sm"
                onClick={() => setBannerDismissed(true)}
              >
                Dismiss
              </button>
            </div>
          )}
          <div className="rg-layout__content">
            <Outlet />
          </div>
        </main>
        <Footer />
      </div>
      {hasJobs && <JobConsole />}
    </div>
  );
}
