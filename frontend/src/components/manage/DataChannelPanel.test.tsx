/**
 * The data channel panel lists channels everywhere and adds, syncs, or removes
 * them only where the instance allows (`data_channels`). Every channel is
 * untrusted for now, so the add form carries a standing warning and each row
 * says so.
 */

import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { describe, expect, it } from 'vitest';
import { DATA_CHANNEL_WARNING, DataChannelPanel } from './DataChannelPanel';
import { useInvalidationBridge } from '../../hooks/useInvalidationBridge';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import { server } from '../../test/server';
import { LOCAL_API } from '../../test/handlers';

function InvalidationBridge() {
  useInvalidationBridge();
  return null;
}

const channel = {
  name: 'registry',
  type: 'https',
  index_address: 'https://example.org/index.yaml',
  description: null,
  credentials_set: false,
  trusted: false,
};

function seedChannels(channels: unknown[]) {
  server.use(
    http.get(`${LOCAL_API}/remote/data_channels`, () => HttpResponse.json({ channels })),
  );
}

describe('DataChannelPanel', () => {
  it('shows the warning and posts the new channel with the action header', async () => {
    let listCalls = 0;
    server.use(
      http.get(`${LOCAL_API}/remote/data_channels`, () => {
        listCalls += 1;
        return HttpResponse.json({ channels: listCalls > 1 ? [channel] : [] });
      }),
    );
    let headerSeen: string | null = null;
    let body: unknown;
    server.use(
      http.post(`${LOCAL_API}/actions/data_channels`, async ({ request }) => {
        headerSeen = request.headers.get('x-refgenie-action');
        body = await request.json();
        return HttpResponse.json({
          ok: true,
          message: `Data channel 'registry' added. ${DATA_CHANNEL_WARNING}`,
          data: { name: 'registry', trusted: false, sync: null },
        });
      }),
    );

    renderWithProviders(
      <>
        <InvalidationBridge />
        <DataChannelPanel />
      </>,
      { route: '/manage' },
    );
    await screen.findByText(/No data channels/);
    expect(screen.getByText(DATA_CHANNEL_WARNING)).toBeInTheDocument();

    await userEvent.type(screen.getByLabelText(/^Name \*/), 'registry');
    await userEvent.type(screen.getByLabelText(/^Index URL \*/), 'https://example.org/index.yaml');
    await userEvent.click(screen.getByRole('button', { name: 'Add channel' }));

    await waitFor(() =>
      expect(body).toEqual({
        name: 'registry',
        index_address: 'https://example.org/index.yaml',
        description: null,
      }),
    );
    expect(headerSeen).not.toBeNull();
    // The list refetches after the add and the new row shows up, unverified.
    await screen.findByRole('cell', { name: 'registry' });
    expect(screen.getByText('not verified')).toBeInTheDocument();
  });

  it('rejects a non-http index before calling the backend', async () => {
    seedChannels([]);
    let posted = false;
    server.use(
      http.post(`${LOCAL_API}/actions/data_channels`, () => {
        posted = true;
        return HttpResponse.json({ ok: true, message: '', data: null });
      }),
    );

    renderWithProviders(<DataChannelPanel />, { route: '/manage' });
    await screen.findByText(/No data channels/);
    await userEvent.type(screen.getByLabelText(/^Name \*/), 'local');
    await userEvent.type(screen.getByLabelText(/^Index URL \*/), '/etc/index.yaml');
    await userEvent.click(screen.getByRole('button', { name: 'Add channel' }));

    expect(await screen.findByRole('alert')).toHaveTextContent(/http\(s\) URL/);
    expect(posted).toBe(false);
  });

  it('syncs a channel from its row', async () => {
    seedChannels([channel]);
    let synced = false;
    server.use(
      http.post(`${LOCAL_API}/actions/data_channels/registry/sync`, () => {
        synced = true;
        return HttpResponse.json({
          ok: true,
          message: "Synced 'registry': 2 new item(s)",
          data: {
            sync: {
              channel: 'registry',
              asset_classes_added: 1,
              asset_classes_skipped: 0,
              asset_classes_failed: 0,
              recipes_added: 1,
              recipes_skipped: 0,
              recipes_failed: 0,
              errors: [],
            },
          },
        });
      }),
    );

    renderWithProviders(<DataChannelPanel />, { route: '/manage' });
    await userEvent.click(await screen.findByRole('button', { name: 'Sync' }));
    await waitFor(() => expect(synced).toBe(true));
  });

  it('removes a channel after confirming', async () => {
    seedChannels([channel]);
    let deleted = false;
    server.use(
      http.delete(`${LOCAL_API}/actions/data_channels/registry`, () => {
        deleted = true;
        return HttpResponse.json({ ok: true, message: 'Removed.', data: null });
      }),
    );

    renderWithProviders(<DataChannelPanel />, { route: '/manage' });
    await userEvent.click(await screen.findByRole('button', { name: 'Remove' }));
    await screen.findByRole('heading', { name: 'Remove data channel' });
    const confirm = screen
      .getAllByRole('button', { name: 'Remove' })
      .find((button) => button.closest('.modal__footer'));
    await userEvent.click(confirm!);
    await waitFor(() => expect(deleted).toBe(true));
  });

  it('lists channels read-only on an instance without data_channels', async () => {
    seedChannels([channel]);
    renderWithProviders(<DataChannelPanel />, { route: '/manage', config: serverConfig });
    await screen.findByRole('cell', { name: 'registry' });
    expect(screen.getByText('not verified')).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'Add channel' })).not.toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'Sync' })).not.toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'Remove' })).not.toBeInTheDocument();
    expect(screen.queryByText(DATA_CHANNEL_WARNING)).not.toBeInTheDocument();
  });
});
