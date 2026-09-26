/**
 * The first-class connect affordance on the landing page while disconnected
 * (bridge plan B4), and the Safari fallback card (B7) on WebKit-only
 * browsers, where the probe is never attempted because WebKit blocks
 * HTTPS-page-to-localhost fetches outright (WebKit bug 171934).
 *
 * Renders nothing once connected, and nothing at all when this page IS the
 * local dash.
 */

import { useState } from 'react';
import { CopyButton } from '../common/CopyButton';
import { ExternalLink } from '../common/ExternalLink';
import { ConnectDialog } from './ConnectDialog';
import { useBridge } from '../../hooks/useBridge';
import {
  dismissSafariCard,
  isSafariCardDismissed,
} from '../../services/bridge/persistence';

export function ConnectCard() {
  const bridge = useBridge();
  const [dialogOpen, setDialogOpen] = useState(false);
  const [safariDismissed, setSafariDismissed] = useState(isSafariCardDismissed());

  if (bridge.status === 'connected' || bridge.status === 'self') return null;

  if (bridge.status === 'blocked') {
    if (safariDismissed) return null;
    const localUrl = `http://localhost:${bridge.port}`;
    return (
      <section className="rg-card">
        <div className="rg-card__body flex items-start justify-between gap-4 flex-wrap">
          <div className="flex flex-col gap-2">
            <h2 className="text-xl font-semibold">Using refgenie locally?</h2>
            <p className="text-sm">
              Safari blocks pages served over HTTPS from talking to a server on your own
              computer (
              <ExternalLink href="https://bugs.webkit.org/show_bug.cgi?id=171934">
                WebKit bug 171934
              </ExternalLink>
              ). Open your local refgenie directly — it runs the same interface:
            </p>
            <span className="flex items-center gap-2">
              <code className="rg-code rg-code--inline">{localUrl}</code>
              <CopyButton value={localUrl} label={`Copy: ${localUrl}`} />
            </span>
          </div>
          <button
            type="button"
            className="rg-btn rg-btn--sm"
            onClick={() => {
              dismissSafariCard();
              setSafariDismissed(true);
            }}
          >
            Don&rsquo;t show again
          </button>
        </div>
      </section>
    );
  }

  return (
    <section className="rg-card">
      <div className="rg-card__body flex items-center justify-between gap-4 flex-wrap">
        <div className="flex flex-col gap-1">
          <h2 className="text-xl font-semibold">Using refgenie locally?</h2>
          <p className="rg-muted text-sm">
            Connect this page to your running{' '}
            <code className="rg-code rg-code--inline">refgenie dash</code> to see which
            assets you already have and browse your local genomes.
          </p>
        </div>
        <button
          type="button"
          className="rg-btn rg-btn--primary"
          onClick={() => setDialogOpen(true)}
        >
          Connect
        </button>
      </div>
      <ConnectDialog open={dialogOpen} onClose={() => setDialogOpen(false)} />
    </section>
  );
}
