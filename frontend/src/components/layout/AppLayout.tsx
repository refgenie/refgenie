import { useState } from 'react';
import { Outlet } from 'react-router-dom';
import { NavBar } from './NavBar';
import { Footer } from './Footer';
import { BridgeControl } from '../bridge/BridgeControl';
import { Icon } from '../common/Icon';
import { PROJECT_LINKS } from '../../services/projectLinks';
import { BridgeJobsDock } from '../bridge/BridgeJobsDock';
import { JobConsole } from '../jobs/JobConsole';
import { cn } from '../../utils/cn';
import { useUiConfig } from '../../hooks/useUiConfig';
import { useCapability } from '../../hooks/useCapability';
import { useJobEvents } from '../../hooks/useJobEvents';
import { useBridgeAutoConnect } from '../../hooks/useBridgeAutoConnect';
import { useBridgeJobEntries } from '../../hooks/useBridgeJobPolling';

export function AppLayout() {
  const config = useUiConfig();
  const hasJobs = useCapability('jobs');
  const [bannerDismissed, setBannerDismissed] = useState(false);
  const isPublic = config.mode !== 'local';
  // The bridge jobs dock needs the same bottom clearance as the job console.
  const hasBridgeJobs = useBridgeJobEntries().length > 0 && isPublic;

  // The one job event stream, mounted once.
  useJobEvents({ enabled: hasJobs });
  // Re-probe a remembered bridge connection (never on a local dash).
  useBridgeAutoConnect(config.mode !== 'local');

  return (
    <div className={cn('rg-layout', (hasJobs || hasBridgeJobs) && 'rg-layout--docked')}>
      <header className="rg-layout__header">
        <div className="rg-layout__bar">
          <NavBar />
          <div className="flex items-center gap-4">
            {config.mode !== 'local' && <BridgeControl />}
            <a
              className="rg-nav__item"
              href={PROJECT_LINKS.docs}
              target="_blank"
              rel="noreferrer"
            >
              <Icon name="book" className="rg-icon--leading" />
              Docs
            </a>
          </div>
        </div>
      </header>
      <div className="rg-layout__body">
        <main className="rg-layout__main">
          {config.degraded && !bannerDismissed && (
            <div
              className="rg-banner rg-banner--warning rg-banner--strip"
              role="status"
            >
              <span>
                Could not read{' '}
                <code className="rg-code rg-code--inline">/service-info</code>, so this
                instance&rsquo;s capabilities are unknown. The interface is read-only.
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
      {isPublic && <BridgeJobsDock />}
    </div>
  );
}
