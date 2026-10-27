/**
 * Data channels: where this instance gets its recipes and asset classes.
 *
 * Every channel is untrusted for now. Nothing verifies who publishes an index,
 * and a recipe synced from one is shell that runs here at build time, so the
 * add form carries a standing warning rather than a checkbox to click past.
 * `trusted` is read off each row so the badge changes by itself when the
 * backend starts vouching for a channel.
 *
 * Add syncs right away by default: a channel that is recorded but never synced
 * is a row with nothing behind it. Sync is offered per row to refresh.
 */

import { useId, useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useDataChannels } from '../../hooks/queries/useRemote';
import { useToast } from '../../hooks/useToast';
import { addDataChannel, removeDataChannel, syncDataChannel } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ApiError } from '../../services/http';
import type { DataChannelSyncReport } from '../../services/contracts';
import { Badge } from '../common/Badge';
import { ConfirmModal } from '../common/ConfirmModal';
import { DataTable } from '../common/DataTable';
import { FormField } from '../common/FormField';
import { EmptyState } from '../common/states';
import type { Column } from '../common/DataTable';
import type { LocalDataChannel } from '../../types/api';

export const DATA_CHANNEL_WARNING =
  'Data channels are not verified. Recipes from a channel run shell commands on this machine when you build. Only add channels you trust.';

function isHttpUrl(value: string): boolean {
  try {
    const parsed = new URL(value);
    return parsed.protocol === 'http:' || parsed.protocol === 'https:';
  } catch {
    return false;
  }
}

function describeSync(report: DataChannelSyncReport | undefined): string {
  if (!report) return '';
  const added = report.asset_classes_added + report.recipes_added;
  const failed = report.asset_classes_failed + report.recipes_failed;
  return failed ? `${added} new item(s), ${failed} failed.` : `${added} new item(s).`;
}

function ChannelActions({ channel }: { channel: LocalDataChannel }) {
  const client = useLocalApiClient();
  const toast = useToast();
  const [confirmRemove, setConfirmRemove] = useState(false);
  const [syncing, setSyncing] = useState(false);

  const sync = async () => {
    setSyncing(true);
    try {
      const result = await syncDataChannel(client, channel.name);
      invalidate(INVALIDATION_FOR_ACTION['data_channel.sync']);
      const report = result.data?.sync as DataChannelSyncReport | undefined;
      const summary = describeSync(report);
      if (report && report.errors.length) toast.error(`${result.message}. ${summary}`);
      else toast.success(`${result.message}. ${summary}`);
    } catch (caught) {
      toast.error(caught instanceof ApiError ? caught.detail : String(caught));
    } finally {
      setSyncing(false);
    }
  };

  return (
    <div className="flex gap-2">
      <button
        type="button"
        className="rg-btn rg-btn--sm"
        disabled={syncing}
        onClick={() => void sync()}
      >
        {syncing ? 'Syncing…' : 'Sync'}
      </button>
      <button
        type="button"
        className="rg-btn rg-btn--sm rg-btn--danger"
        onClick={() => setConfirmRemove(true)}
      >
        Remove
      </button>

      <ConfirmModal
        isOpen={confirmRemove}
        onClose={() => setConfirmRemove(false)}
        title="Remove data channel"
        destructive
        confirmLabel="Remove"
        body={
          <p>
            Forget <code className="rg-code rg-code--inline">{channel.name}</code>? Recipes and
            asset classes already synced from it stay registered.
          </p>
        }
        onConfirm={async () => {
          await removeDataChannel(client, channel.name);
          invalidate(INVALIDATION_FOR_ACTION['data_channel.remove']);
          toast.success(`Data channel ${channel.name} removed.`);
        }}
      />
    </div>
  );
}

function AddChannelForm({ knownNames }: { knownNames: readonly string[] }) {
  const client = useLocalApiClient();
  const toast = useToast();
  const nameId = useId();
  const urlId = useId();
  const descriptionId = useId();

  const [name, setName] = useState('');
  const [url, setUrl] = useState('');
  const [description, setDescription] = useState('');
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const submit = async (event: React.FormEvent) => {
    event.preventDefault();
    const channelName = name.trim();
    const indexAddress = url.trim();
    if (knownNames.includes(channelName)) {
      setError(`A data channel named ${channelName} already exists.`);
      return;
    }
    if (!isHttpUrl(indexAddress)) {
      setError('Enter a full http(s) URL to the channel’s index.yaml.');
      return;
    }
    setBusy(true);
    setError(null);
    try {
      const result = await addDataChannel(client, {
        name: channelName,
        index_address: indexAddress,
        description: description.trim() || null,
      });
      invalidate(INVALIDATION_FOR_ACTION['data_channel.add']);
      const report = result.data?.sync as DataChannelSyncReport | undefined;
      if (report && report.errors.length) {
        toast.error(`Added ${channelName}, but the sync had failures. ${describeSync(report)}`);
      } else {
        toast.success(`Data channel ${channelName} added. ${describeSync(report)}`.trim());
      }
      setName('');
      setUrl('');
      setDescription('');
    } catch (caught) {
      setError(caught instanceof ApiError ? caught.detail : String(caught));
    } finally {
      setBusy(false);
    }
  };

  return (
    <form className="flex flex-col gap-3" noValidate onSubmit={(event) => void submit(event)}>
      <h3 className="text-lg font-semibold">Add a data channel</h3>

      <FormField htmlFor={nameId} label="Name" required hint="A short label, for example refgenie.">
        <input
          id={nameId}
          className="rg-field__input"
          type="text"
          value={name}
          onChange={(event) => {
            setName(event.target.value);
            setError(null);
          }}
        />
      </FormField>

      <FormField
        htmlFor={urlId}
        label="Index URL"
        required
        error={error ?? undefined}
        hint="For example https://refgenie.github.io/refgenie-registry/index.yaml"
      >
        <input
          id={urlId}
          className="rg-field__input rg-field__input--mono"
          type="url"
          value={url}
          placeholder="https://…/index.yaml"
          onChange={(event) => {
            setUrl(event.target.value);
            setError(null);
          }}
        />
      </FormField>

      <FormField htmlFor={descriptionId} label="Description">
        <input
          id={descriptionId}
          className="rg-field__input"
          type="text"
          value={description}
          onChange={(event) => setDescription(event.target.value)}
        />
      </FormField>

      <div className="rg-banner rg-banner--warning" role="status">
        <span>{DATA_CHANNEL_WARNING}</span>
      </div>

      <div className="rg-form-actions">
        <button
          type="submit"
          className="rg-btn rg-btn--primary"
          disabled={busy || !name.trim() || !url.trim()}
        >
          {busy ? 'Adding…' : 'Add channel'}
        </button>
      </div>
    </form>
  );
}

export function DataChannelPanel() {
  const canWrite = useCapability('data_channels');
  const channels = useDataChannels();

  const rows = channels.data?.channels ?? [];

  const columns: Array<Column<LocalDataChannel>> = [
    {
      key: 'name',
      header: 'Channel',
      render: (channel) => (
        <div className="flex flex-col">
          <span className="font-medium">{channel.name}</span>
          {channel.description && <span className="rg-muted text-sm">{channel.description}</span>}
        </div>
      ),
    },
    {
      key: 'index_address',
      header: 'Index',
      render: (channel) => (
        <code className="rg-code rg-code--inline">{channel.index_address}</code>
      ),
    },
    { key: 'type', header: 'Type', render: (channel) => channel.type },
    {
      key: 'trusted',
      header: 'Trust',
      render: (channel) =>
        channel.trusted ? (
          <Badge variant="local">verified</Badge>
        ) : (
          <Badge title={DATA_CHANNEL_WARNING}>not verified</Badge>
        ),
    },
    ...(canWrite
      ? [
          {
            key: 'actions',
            header: '',
            render: (channel: LocalDataChannel) => <ChannelActions channel={channel} />,
          },
        ]
      : []),
  ];

  return (
    <section className="flex flex-col gap-4">
      <h2 className="text-xl font-semibold">Data channels</h2>

      <DataTable
        caption="Data channels"
        columns={columns}
        rows={rows}
        rowKey={(channel) => channel.name}
        loading={channels.isPending}
        error={channels.error}
        onRetry={() => channels.refetch()}
        empty={
          <EmptyState message="No data channels. Add one to get recipes and asset classes to build from." />
        }
      />

      {canWrite && <AddChannelForm knownNames={rows.map((channel) => channel.name)} />}
    </section>
  );
}
