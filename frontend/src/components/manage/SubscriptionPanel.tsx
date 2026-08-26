/**
 * Server subscriptions.
 *
 * The list and the reachability badge both come from `/v1/remote/servers`,
 * which the browse UI already probes — probing again from here would double the
 * network cost and could disagree with what the remote catalogue is showing.
 *
 * Subscribe and unsubscribe invalidate `servers` AND `remote`, which forces the
 * remote catalogue to refetch; otherwise it keeps listing assets from a server
 * the user just dropped.
 */

import { useState } from 'react';
import { useRemoteServers } from '../../hooks/queries/useRemote';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { subscribe, unsubscribe } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ApiError } from '../../services/http';
import { DataTable } from '../common/DataTable';
import { Badge } from '../common/Badge';
import { ConfirmModal } from '../common/ConfirmModal';
import { FormField } from '../common/FormField';
import { EmptyState } from '../common/states';
import type { Column } from '../common/DataTable';
import type { RemoteServer } from '../../types/api';

function isValidUrl(value: string): boolean {
  try {
    const parsed = new URL(value);
    return parsed.protocol === 'http:' || parsed.protocol === 'https:';
  } catch {
    return false;
  }
}

export function SubscriptionPanel() {
  const canWrite = useCapability('subscriptions');
  const client = useLocalApiClient();
  const toast = useToast();
  const servers = useRemoteServers();

  const [url, setUrl] = useState('');
  const [reset, setReset] = useState(false);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [pending, setPending] = useState<RemoteServer | null>(null);

  const rows = (servers.data?.servers ?? []).filter((server) => server.subscribed);

  const submit = async (event: React.FormEvent) => {
    event.preventDefault();
    if (!isValidUrl(url)) {
      setError('Enter a full http(s) URL.');
      return;
    }
    setBusy(true);
    setError(null);
    try {
      // Always the list form, and never empty: the model sets min_length=1.
      await subscribe(client, { server_urls: [url.trim()], reset });
      invalidate(INVALIDATION_FOR_ACTION.subscribe);
      toast.success(`Subscribed to ${url.trim()}.`);
      setUrl('');
      setReset(false);
    } catch (caught) {
      setError(caught instanceof ApiError ? caught.detail : String(caught));
    } finally {
      setBusy(false);
    }
  };

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
            render: (server: RemoteServer) => (
              <button
                type="button"
                className="rg-btn rg-btn--sm rg-btn--danger"
                onClick={() => setPending(server)}
              >
                Unsubscribe
              </button>
            ),
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

      {canWrite && (
        <form className="flex flex-col gap-3" onSubmit={(event) => void submit(event)}>
          <FormField
            htmlFor="subscribe-url"
            label="Subscribe to a server"
            error={error ?? undefined}
            hint="For example https://refgenomes.databio.org"
          >
            <input
              id="subscribe-url"
              className="rg-field__input rg-field__input--mono"
              type="url"
              value={url}
              placeholder="https://refgenomes.databio.org"
              onChange={(event) => setUrl(event.target.value)}
            />
          </FormField>

          <label className="rg-check" htmlFor="subscribe-reset">
            <input
              id="subscribe-reset"
              type="checkbox"
              checked={reset}
              onChange={(event) => setReset(event.target.checked)}
            />
            <span>and replace existing subscriptions</span>
          </label>

          <div className="rg-form-actions">
            <button type="submit" className="rg-btn rg-btn--primary" disabled={busy || !url}>
              {busy ? 'Subscribing…' : 'Subscribe'}
            </button>
          </div>
        </form>
      )}

      <ConfirmModal
        isOpen={pending !== null}
        onClose={() => setPending(null)}
        title="Unsubscribe"
        destructive
        confirmLabel="Unsubscribe"
        body={
          <p>
            Stop using <code className="rg-code rg-code--inline">{pending?.url}</code>? Assets
            already pulled from it stay; nothing new can be pulled from it.
          </p>
        }
        onConfirm={async () => {
          if (!pending) return;
          await unsubscribe(client, { server_urls: [pending.url] });
          invalidate(INVALIDATION_FOR_ACTION.unsubscribe);
          toast.success(`Unsubscribed from ${pending.url}.`);
          setPending(null);
        }}
      />
    </section>
  );
}
