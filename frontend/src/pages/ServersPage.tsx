/**
 * The refgenie servers this install pulls from: the subscription list with
 * reachability, plus subscribe and unsubscribe on a local dash.
 *
 * `?subscribe=<url>` opens with the form filled in (never submitted). Failed
 * pulls link here with the server they targeted, so "not subscribed" is one
 * click from fixed.
 *
 * Under the bridge's `/local` scope the page is read-only: subscribing changes
 * the local config, which only the local dash itself may do. The page links
 * there instead of offering a form that the local refgenie would refuse.
 */

import { useSearchParams } from 'react-router-dom';
import { useBridge } from '../hooks/useBridge';
import { useCapability } from '../hooks/useCapability';
import { useRouteBase } from '../hooks/useRouteBase';
import { MiniHero } from '../components/layout/MiniHero';
import { ExternalLink } from '../components/common/ExternalLink';
import { SubscriptionPanel } from '../components/servers/SubscriptionPanel';
import { SUBSCRIBE_PARAM, localServersUrl } from '../components/servers/links';
import { NotAvailablePage } from './NotAvailablePage';

export function ServersPage() {
  const canBrowse = useCapability('remote_browse');
  const canSubscribe = useCapability('subscriptions');
  const routeBase = useRouteBase();
  const bridge = useBridge();
  const [params] = useSearchParams();
  const prefill = params.get(SUBSCRIBE_PARAM)?.trim() || undefined;
  const bridged = routeBase === '/local';

  // The list itself is `/v1/remote/servers`, local mode only. Inside `/local`
  // the masked config hides `remote_browse`, but the bridged local refgenie
  // serves that read, so the scope alone is enough.
  if (!canBrowse && !canSubscribe && !bridged) {
    return (
      <NotAvailablePage
        title="Servers"
        reason="This instance does not pull from other servers."
      />
    );
  }

  return (
    <div className="flex flex-col gap-6">
      <MiniHero
        title="Servers"
        lede={
          <>
            A server is a refgenie that hosts pre-built assets for download. Your refgenie pulls
            only from the servers it subscribes to; a new install subscribes to{' '}
            <code className="rg-code rg-code--inline">https://api.refgenie.org</code>.
          </>
        }
      />

      <div className="flex flex-col gap-6 max-w-screen-lg mx-auto w-full">
        {bridged && bridge.baseUrl && (
          <div className="rg-banner rg-banner--info" role="note">
            <span>
              This page cannot change your subscriptions. To add or remove a server,{' '}
              <ExternalLink href={localServersUrl(bridge.baseUrl, prefill)}>
                open Servers on your local refgenie
              </ExternalLink>
              .
            </span>
          </div>
        )}
        <SubscriptionPanel prefillUrl={prefill} />
      </div>
    </div>
  );
}
