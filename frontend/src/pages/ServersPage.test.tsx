/**
 * The Servers page: the subscription list on a local dash, the `?subscribe=`
 * prefill that failed pulls link to, and the read-only `/local` view on the
 * public site, which links to the local dash instead of subscribing.
 */

import { beforeEach, describe, expect, it } from 'vitest';
import { screen, within } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { Route, Routes } from 'react-router-dom';
import { ServersPage } from './ServersPage';
import { LocalScope } from '../components/bridge/LocalScope';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';
import { server } from '../test/server';
import { remoteServersFixture } from '../test/fixtures';
import { useBridgeStore } from '../stores/bridgeStore';
import { coerceCapabilities } from '../services/capabilities';
import { CAPABILITY_KEYS } from '../types/ui';
import { serversPath } from '../components/servers/links';
import type { PingResponse } from '../services/bridge/contract';

const BRIDGE_ORIGIN = 'http://localhost:8080';

const ping: PingResponse = {
  service: 'refgenie',
  bridge_version: 1,
  mode: 'local',
  refgenie_version: '1.0.0',
  instance_id: 'instance-1',
  instance_label: 'my laptop refgenie',
  bridge_mode: 'full',
  action_header: 'X-Refgenie-Action',
  capabilities: coerceCapabilities(Object.fromEntries(CAPABILITY_KEYS.map((key) => [key, true]))),
  bridge: { actions_cross_origin: true },
};

beforeEach(() => {
  localStorage.clear();
  useBridgeStore.getState().reset();
});

describe('ServersPage (local dash)', () => {
  it('lists subscriptions with reachability, unsubscribe, and an empty form', async () => {
    renderWithProviders(<ServersPage />, { route: '/servers' });

    expect(screen.getByRole('heading', { level: 1, name: 'Servers' })).toBeInTheDocument();
    const table = await screen.findByRole('table', { name: 'Subscribed servers' });
    expect(within(table).getByText('https://api.refgenie.org')).toBeInTheDocument();
    expect(within(table).getByText('reachable')).toBeInTheDocument();
    expect(within(table).getByText('unreachable')).toBeInTheDocument();
    expect(within(table).getAllByRole('button', { name: 'Unsubscribe' })).toHaveLength(2);
    expect(screen.getByLabelText('Subscribe to a server')).toHaveValue('');
  });

  it('opens the form prefilled from ?subscribe=, without submitting it', async () => {
    let posted = false;
    server.use(
      http.post('/v1/actions/subscriptions', () => {
        posted = true;
        return HttpResponse.json({ ok: true, message: '', data: null });
      }),
    );
    renderWithProviders(<ServersPage />, { route: serversPath('https://example.org') });

    const field = screen.getByLabelText('Subscribe to a server');
    expect(field).toHaveValue('https://example.org');
    expect(field).toHaveFocus();
    await screen.findByRole('table', { name: 'Subscribed servers' });
    expect(posted).toBe(false);
  });

  it('says so when nothing is subscribed', async () => {
    server.use(http.get('/v1/remote/servers', () => HttpResponse.json({ servers: [] })));
    renderWithProviders(<ServersPage />, { route: '/servers' });
    expect(await screen.findByText(/no servers subscribed/i)).toBeInTheDocument();
  });
});

describe('ServersPage under /local (public site, bridged)', () => {
  it('lists the local subscriptions read-only and links to the local dash to subscribe', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    server.use(
      http.get(`${BRIDGE_ORIGIN}/v1/remote/servers`, () =>
        HttpResponse.json(remoteServersFixture),
      ),
    );
    renderWithProviders(
      <Routes>
        <Route
          path="/local/servers"
          element={
            <LocalScope>
              <ServersPage />
            </LocalScope>
          }
        />
      </Routes>,
      { config: serverConfig, route: '/local/servers?subscribe=https%3A%2F%2Fexample.org' },
    );

    const table = await screen.findByRole('table', { name: 'Subscribed servers' });
    expect(within(table).getByText('https://api.refgenie.org')).toBeInTheDocument();
    // No write controls: subscribing is never a bridge action.
    expect(screen.queryByRole('button', { name: 'Unsubscribe' })).toBeNull();
    expect(screen.queryByRole('button', { name: 'Subscribe' })).toBeNull();
    const link = screen.getByRole('link', { name: /open servers on your local refgenie/i });
    expect(link).toHaveAttribute(
      'href',
      `${BRIDGE_ORIGIN}/servers?subscribe=${encodeURIComponent('https://example.org')}`,
    );
    // The scope's own nav marks Servers as current.
    const nav = screen.getByRole('navigation', { name: 'Local refgenie' });
    expect(within(nav).getByRole('link', { name: 'Servers' })).toHaveClass(
      'rg-nav__item--active',
    );
  });

  it('highlights Genomes on bare /local', () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    renderWithProviders(
      <LocalScope>
        <span />
      </LocalScope>,
      { config: serverConfig, route: '/local' },
    );
    const nav = screen.getByRole('navigation', { name: 'Local refgenie' });
    expect(within(nav).getByRole('link', { name: 'Genomes' })).toHaveClass(
      'rg-nav__item--active',
    );
    expect(within(nav).getByRole('link', { name: 'Species' })).not.toHaveClass(
      'rg-nav__item--active',
    );
  });
});
