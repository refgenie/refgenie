/**
 * Top-bar bridge control: a Connect button while disconnected, a status pill
 * with local-browse and disconnect once connected. Hidden on the local dash
 * itself and reduced to one line on WebKit. Kept on one line, so counts go in the tooltip.
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
          Not available in Safari. Instead,{' '}
          <ExternalLink href={`http://localhost:${bridge.port}`}>
            open local refgenie
          </ExternalLink>
          .
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
