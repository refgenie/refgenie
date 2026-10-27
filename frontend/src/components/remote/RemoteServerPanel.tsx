import { useCapability } from '../../hooks/useCapability';
import { SubscribeForm } from '../actions/SubscribeForm';
import { UnsubscribeButton } from '../actions/UnsubscribeButton';
import { cn } from '../../utils/cn';
import type { RemoteServer } from '../../types/api';

export interface RemoteServerPanelProps {
  servers: RemoteServer[] | undefined;
  selected: string | undefined;
  onSelect: (url: string) => void;
  /** A server was just added, for a caller that wants to browse it at once. */
  onSubscribed?: (url: string) => void;
  /** A server is gone, for a caller holding it as its selection. */
  onUnsubscribed?: (url: string) => void;
}

/**
 * Server picker, and the place a subscription is added or dropped.
 *
 * Per-server reachability and error text are rendered inline, never a silently
 * empty list. An instance with no subscriptions is
 * where the Remote page starts its life, so the form is right here rather than
 * a sentence naming a CLI command; when the instance cannot subscribe over
 * HTTP the CLI command is still the honest answer.
 */
export function RemoteServerPanel({
  servers,
  selected,
  onSelect,
  onSubscribed,
  onUnsubscribed,
}: RemoteServerPanelProps) {
  const canSubscribe = useCapability('subscriptions');
  const rows = servers ?? [];

  return (
    <div className="flex flex-col gap-4">
      {rows.length === 0 ? (
        <p className="rg-muted text-sm">
          {canSubscribe ? (
            'No servers subscribed. Add one below to browse and pull what it holds.'
          ) : (
            <>
              No servers configured. Add one with{' '}
              <code className="rg-code rg-code--inline">refgenie subscribe</code>.
            </>
          )}
        </p>
      ) : (
        <ul className="flex flex-col gap-2">
          {rows.map((server) => (
            <li key={server.url}>
              <div className="flex items-center gap-2">
                <button
                  type="button"
                  className={cn('rg-btn flex-1', server.url === selected && 'rg-btn--primary')}
                  onClick={() => onSelect(server.url)}
                  aria-current={server.url === selected ? 'true' : undefined}
                >
                  <span className="flex-1 text-left">{server.url}</span>
                  <span className="text-xs">
                    {server.reachable ? 'reachable' : 'unreachable'}
                    {server.subscribed ? ' · subscribed' : ''}
                  </span>
                </button>
                <UnsubscribeButton url={server.url} onUnsubscribed={onUnsubscribed} />
              </div>
              {server.error && (
                <p className="rg-banner rg-banner--warning mt-1 text-xs" role="alert">
                  {server.error}
                </p>
              )}
            </li>
          ))}
        </ul>
      )}

      <SubscribeForm
        knownUrls={rows.map((server) => server.url)}
        label={rows.length === 0 ? 'Subscribe to a server' : 'Subscribe to another server'}
        onSubscribed={onSubscribed}
      />
    </div>
  );
}
