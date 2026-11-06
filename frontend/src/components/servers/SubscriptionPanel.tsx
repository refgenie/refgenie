/**
 * Server subscriptions: the body of the Servers page. The list and
 * reachability come from `/v1/remote/servers`. Add and remove hide themselves
 * without the `subscriptions` capability.
 */

import { useRemoteServers } from '../../hooks/queries/useRemote';
import { useCapability } from '../../hooks/useCapability';
import { SubscribeForm } from '../actions/SubscribeForm';
import { UnsubscribeButton } from '../actions/UnsubscribeButton';
import { DataTable } from '../common/DataTable';
import { Badge } from '../common/Badge';
import { EmptyState } from '../common/states';
import type { Column } from '../common/DataTable';
import type { RemoteServer } from '../../types/api';

export interface SubscriptionPanelProps {
  /** Prefill the subscribe form, from the page's `?subscribe=` parameter. */
  prefillUrl?: string;
}

export function SubscriptionPanel({ prefillUrl }: SubscriptionPanelProps) {
  const canWrite = useCapability('subscriptions');
  const servers = useRemoteServers();

  const rows = (servers.data?.servers ?? []).filter((server) => server.subscribed);

  const columns: Array<Column<RemoteServer>> = [
    {
      key: 'url',
      header: 'Server',
      render: (server) => <code className="rg-code rg-code--inline">{server.url}</code>,
    },
    {
      key: 'reachable',
      header: 'Status',
      render: (server) =>
        server.reachable ? (
          <Badge variant="local">reachable</Badge>
        ) : (
          <span className="rg-muted" title={server.error ?? undefined}>
            unreachable
          </span>
        ),
    },
    ...(canWrite
      ? [
          {
            key: 'actions',
            header: '',
            render: (server: RemoteServer) => <UnsubscribeButton url={server.url} />,
          },
        ]
      : []),
  ];

  return (
    <section className="flex flex-col gap-4">
      <h2 className="text-xl font-semibold">Subscribed servers</h2>

      <DataTable
        caption="Subscribed servers"
        columns={columns}
        rows={rows}
        rowKey={(server) => server.url}
        loading={servers.isPending}
        error={servers.error}
        onRetry={() => servers.refetch()}
        empty={<EmptyState message="No servers subscribed. Remote browse and pull need one." />}
      />

      <SubscribeForm
        // Keyed on the prefill so a new `?subscribe=` link resets the field.
        key={prefillUrl ?? ''}
        knownUrls={rows.map((server) => server.url)}
        allowReset
        initialUrl={prefillUrl}
        autoFocus={!!prefillUrl}
      />
    </section>
  );
}
