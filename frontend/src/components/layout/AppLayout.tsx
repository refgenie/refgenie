import { useState } from 'react';
import { Outlet } from 'react-router-dom';
import { NavBar } from './NavBar';
import { Footer } from './Footer';
import { JobConsole } from '../jobs/JobConsole';
import { cn } from '../../utils/cn';
import { useUiConfig } from '../../hooks/useUiConfig';
import { useCapability } from '../../hooks/useCapability';
import { useInvalidationBridge } from '../../hooks/useInvalidationBridge';
import { useJobEvents } from '../../hooks/useJobEvents';

export function AppLayout() {
  const config = useUiConfig();
  const hasJobs = useCapability('jobs');
  const [bannerDismissed, setBannerDismissed] = useState(false);

  // The single event transport and the single query-invalidation listener are
  // both mounted here, once. Nothing else in the tree opens a connection.
  useInvalidationBridge();
  useJobEvents({ enabled: hasJobs });

  return (
    <div className={cn('rg-layout', hasJobs && 'rg-layout--docked')}>
      <div className="rg-layout__sidebar">
        <NavBar />
      </div>
      <div className="rg-layout__body">
        <main className="rg-layout__main">
          <div className="rg-layout__content">
            {config.degraded && !bannerDismissed && (
              <div className="rg-banner rg-banner--warning mb-6" role="status">
                <span>
                  Could not read <code className="rg-code rg-code--inline">/service-info</code>,
                  so this instance&rsquo;s capabilities are unknown. The interface is read-only.
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
            <Outlet />
          </div>
        </main>
        <Footer />
      </div>
      {hasJobs && <JobConsole />}
    </div>
  );
}
