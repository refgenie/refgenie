import { cn } from '../../utils/cn';
import type { RemoteServer } from '../../types/api';

export interface RemoteServerPanelProps {
  servers: RemoteServer[] | undefined;
  selected: string | undefined;
  onSelect: (url: string) => void;
}

/**
 * Server picker. Per-server reachability and error text are rendered inline:
 * a silently empty list was the old dash's failure mode.
 */
export function RemoteServerPanel({ servers, selected, onSelect }: RemoteServerPanelProps) {
  if (!servers || servers.length === 0) {
    return (
      <p className="rg-muted text-sm">
        No servers configured. Add one with <code className="rg-code rg-code--inline">refgenie subscribe</code>.
      </p>
    );
  }

  return (
    <ul className="flex flex-col gap-2">
      {servers.map((server) => (
        <li key={server.url}>
          <button
            type="button"
            className={cn('rg-btn w-full', server.url === selected && 'rg-btn--primary')}
            onClick={() => onSelect(server.url)}
            aria-current={server.url === selected ? 'true' : undefined}
          >
            <span className="flex-1 text-left">{server.url}</span>
            <span className="text-xs">
              {server.reachable ? 'reachable' : 'unreachable'}
              {server.subscribed ? ' · subscribed' : ''}
            </span>
          </button>
          {server.error && (
            <p className="rg-banner rg-banner--warning mt-1 text-xs" role="alert">
              {server.error}
            </p>
          )}
        </li>
      ))}
    </ul>
  );
}
