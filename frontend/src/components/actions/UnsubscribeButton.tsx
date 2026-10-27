/**
 * Drop a server subscription.
 *
 * Paired with `SubscribeForm` and used from the same two surfaces. Nothing on
 * disk is touched: assets already pulled from the server stay, only the ability
 * to reach it for new pulls goes away — which is what the confirmation says.
 *
 * It invalidates `servers` AND `remote`, which forces the remote catalogue to
 * refetch; otherwise it keeps listing assets from a server just dropped.
 */

import { useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { unsubscribe } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ConfirmModal } from '../common/ConfirmModal';

export interface UnsubscribeButtonProps {
  url: string;
  /** The URL that is now gone, for a caller holding it as a selection. */
  onUnsubscribed?: (url: string) => void;
}

export function UnsubscribeButton({ url, onUnsubscribed }: UnsubscribeButtonProps) {
  const canSubscribe = useCapability('subscriptions');
  const client = useLocalApiClient();
  const toast = useToast();
  const [open, setOpen] = useState(false);

  if (!canSubscribe) return null;

  return (
    <>
      <button
        type="button"
        className="rg-btn rg-btn--sm rg-btn--danger"
        onClick={() => setOpen(true)}
      >
        Unsubscribe
      </button>

      <ConfirmModal
        isOpen={open}
        onClose={() => setOpen(false)}
        title="Unsubscribe"
        destructive
        confirmLabel="Unsubscribe"
        body={
          <p>
            Stop using <code className="rg-code rg-code--inline">{url}</code>? Assets already
            pulled from it stay; nothing new can be pulled from it.
          </p>
        }
        onConfirm={async () => {
          await unsubscribe(client, { server_urls: [url] });
          invalidate(INVALIDATION_FOR_ACTION.unsubscribe);
          toast.success(`Unsubscribed from ${url}.`);
          onUnsubscribed?.(url);
        }}
      />
    </>
  );
}
