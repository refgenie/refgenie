/**
 * The top bar's bridge surface: a connect control while disconnected, a
 * connection pill (with local-browse link and disconnect) once connected.
 *
 * Suppressed entirely when the page itself is the local dash (`self`, D7) and
 * reduced to one line on WebKit (`blocked` — the landing page's card carries
 * the full explanation there).
 *
 * It is NOT a `navigation.ts` entry: the surface is stateful (connect link →
 * status pill + browse link + disconnect), and a `{to, label}` record cannot
 * express that. It mounts as a sibling of `<NavBar />` inside the bar, pushed
 * to the trailing edge. Everything sits on ONE line: the bar sets the height
 * of every page, so the genome and asset counts ride the badge's tooltip
 * rather than taking a second row.
 */

import { useState } from 'react';
import { Link } from 'react-router-dom';
import { Badge } from '../common/Badge';
import { ExternalLink } from '../common/ExternalLink';
import { ConnectDialog } from './ConnectDialog';
import { useBridge } from '../../hooks/useBridge';

export function BridgeControl() {
  const bridge = useBridge();
  const [dialogOpen, setDialogOpen] = useState(false);

  if (bridge.status === 'self') return null;

  return (
    <div className="flex flex-wrap items-center gap-2 text-sm">
      {bridge.status === 'connected' ? (
        <>
          <Badge
            variant="local"
            title={
              bridge.ping?.counts
                ? `${bridge.ping.counts.genomes} genomes, ${bridge.ping.counts.assets} assets`
                : undefined
            }
          >
            connected
          </Badge>
          <span className="rg-muted">{bridge.ping?.instance_label}</span>
          <Link className="rg-link" to="/local/genomes">
            Browse local
          </Link>
          <button
            type="button"
            className="rg-btn rg-btn--sm rg-btn--bare"
            onClick={() => bridge.disconnect()}
          >
            Disconnect
          </button>
        </>
      ) : bridge.status === 'blocked' ? (
        <p className="rg-muted">
          Not available in Safari —{' '}
          <ExternalLink href={`http://localhost:${bridge.port}`}>
            open local refgenie
          </ExternalLink>
        </p>
      ) : (
        <button
          type="button"
          className="rg-btn rg-btn--sm rg-btn--bare"
          onClick={() => setDialogOpen(true)}
        >
          Connect local refgenie
        </button>
      )}

      <ConnectDialog open={dialogOpen} onClose={() => setDialogOpen(false)} />
    </div>
  );
}
