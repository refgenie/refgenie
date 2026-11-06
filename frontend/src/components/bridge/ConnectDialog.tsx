/**
 * Dialog for connecting to a local refgenie. The probe runs only on the Connect
 * click, so Chrome's Local Network Access prompt appears when the user expects it.
 * The troubleshooting list covers every cause, because the browser reports
 * mixed content, LNA denial, CORS and connection-refused as the same `TypeError`.
 */

import { useEffect, useId, useState } from 'react';
import { BaseModal } from '../common/BaseModal';
import { FormField } from '../common/FormField';
import { useBridge } from '../../hooks/useBridge';

export interface ConnectDialogProps {
  open: boolean;
  onClose: () => void;
}

export function ConnectDialog({ open, onClose }: ConnectDialogProps) {
  const bridge = useBridge();
  const [portText, setPortText] = useState(String(bridge.port));
  const portId = useId();

  // Close automatically once connected.
  useEffect(() => {
    if (open && bridge.status === 'connected') onClose();
  }, [open, bridge.status, onClose]);

  if (!open) return null;

  const parsedPort = Number(portText);
  const portValid =
    Number.isInteger(parsedPort) && parsedPort > 0 && parsedPort < 65536;

  return (
    <BaseModal
      isOpen={open}
      onClose={onClose}
      title="Connect to my local refgenie"
      size="md"
    >
      <div className="flex flex-col gap-4 text-sm">
        <p>
          If you run refgenie on this computer, this page can show your own genomes next
          to the ones on this server, and mark the assets you already have.
        </p>

        <ol className="flex flex-col gap-2 pl-4 list-decimal">
          <li>
            Start the local dashboard:{' '}
            <code className="rg-code rg-code--inline">refgenie dash</code>. It listens
            on port 8080 unless you pass{' '}
            <code className="rg-code rg-code--inline">--port</code>.
          </li>
          <li>
            Click Connect. Your browser may ask to allow access to devices on your local
            network; allow it.
          </li>
          <li>
            Use <strong>Browse local</strong> at the top of the page to see your
            genomes. To pull assets from this page into your local refgenie, start the
            dashboard with{' '}
            <code className="rg-code rg-code--inline">refgenie dash --bridge full</code>
            .
          </li>
        </ol>

        <p className="rg-muted">
          Your browser talks to your local refgenie directly. Nothing about your local
          genomes is sent to this server. Safari does not allow this connection; open{' '}
          <code className="rg-code rg-code--inline">http://localhost:8080</code> there
          instead.
        </p>

        <FormField
          htmlFor={portId}
          label="Port"
          error={portValid ? undefined : 'Enter a port between 1 and 65535.'}
        >
          <input
            id={portId}
            className="rg-field__input"
            type="text"
            inputMode="numeric"
            value={portText}
            onChange={(event) => setPortText(event.target.value)}
          />
        </FormField>

        {bridge.status === 'unsupported' && (
          <div className="rg-banner rg-banner--warning" role="status">
            <span>
              {bridge.unsupportedReason === 'unsupported-version' &&
                'Your local refgenie speaks a newer bridge protocol than this page expects. Open it directly instead.'}
              {bridge.unsupportedReason === 'wrong-service' &&
                'Something answered on that port, but it does not look like refgenie.'}
            </span>
          </div>
        )}

        {bridge.status === 'absent' && bridge.showTroubleshooting && (
          <div className="rg-banner rg-banner--warning" role="status">
            <div className="flex flex-col gap-2">
              <p className="font-semibold">Could not reach a local refgenie. Check:</p>
              <ul className="flex flex-col gap-1 pl-4">
                <li>
                  Is it running? Start it with{' '}
                  <code className="rg-code rg-code--inline">refgenie dash</code>.
                </li>
                <li>
                  Did your browser ask to allow local network access? A denied prompt
                  can be re-granted from the site settings icon in the address bar.
                </li>
                <li>Is it on a different port? Adjust the port field above.</li>
                <li>
                  Is the bridge disabled? Restart with{' '}
                  <code className="rg-code rg-code--inline">
                    refgenie dash --bridge read
                  </code>{' '}
                  (or check{' '}
                  <code className="rg-code rg-code--inline">REFGENIE_BRIDGE_MODE</code>
                  ).
                </li>
              </ul>
            </div>
          </div>
        )}
      </div>

      <BaseModal.Footer>
        <button type="button" className="rg-btn" onClick={onClose}>
          Cancel
        </button>
        <button
          type="button"
          className="rg-btn rg-btn--primary"
          disabled={!portValid || bridge.status === 'connecting'}
          onClick={() => void bridge.connect(parsedPort)}
        >
          {bridge.status === 'connecting' ? 'Connecting…' : 'Connect'}
        </button>
      </BaseModal.Footer>
    </BaseModal>
  );
}
