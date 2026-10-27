/**
 * Server subscriptions.
 *
 * The list and the reachability badge both come from `/v1/remote/servers`,
 * which the browse UI already probes — probing again from here would double the
 * network cost and could disagree with what the remote catalogue is showing.
 *
 * The add and remove affordances themselves live in `components/actions`,
 * because the Remote page's server picker offers the same two operations right
 * next to the list it is already showing.
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

export function SubscriptionPanel() {
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

      <SubscribeForm knownUrls={rows.map((server) => server.url)} allowReset />
    </section>
  );
}
