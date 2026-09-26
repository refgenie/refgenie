import { beforeEach, describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { LocalScope } from './LocalScope';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import { useBridgeStore } from '../../stores/bridgeStore';
import { useApiClient } from '../../hooks/useApiClient';
import { useUiConfig } from '../../hooks/useUiConfig';
import { useRouteBase } from '../../hooks/useRouteBase';
import { coerceCapabilities } from '../../services/capabilities';
import { CAPABILITY_KEYS } from '../../types/ui';
import type { PingResponse } from '../../services/bridge/contract';
import type { CapabilityKey } from '../../types/ui';

/** Everything a `/local/*` page could use to reach the wrong backend. */
function Probe() {
  const config = useUiConfig();
  const client = useApiClient();
  const routeBase = useRouteBase();
  const on = CAPABILITY_KEYS.filter((key) => config.capabilities[key]);
  return (
    <div>
      <span data-testid="base">{client.baseUrl}</span>
      <span data-testid="route-base">{routeBase}</span>
      <span data-testid="caps">{on.join(',')}</span>
      <span data-testid="service">{config.service_name}</span>
    </div>
  );
}

/** A local instance that claims it can do everything. */
const permissivePing: PingResponse = {
  service: 'refgenie',
  bridge_version: 1,
  mode: 'local',
  refgenie_version: '1.0.0',
  instance_id: 'instance-1',
  instance_label: 'my laptop refgenie',
  bridge_mode: 'full',
  action_header: 'X-Refgenie-Action',
  capabilities: coerceCapabilities(
    Object.fromEntries(CAPABILITY_KEYS.map((key) => [key, true])),
  ),
  bridge: { actions_cross_origin: true },
};

/** Read-only surfaces the scope is allowed to keep. */
const ALLOWED: CapabilityKey[] = ['downloads', 'archives', 'seqcol', 'drs'];

beforeEach(() => {
  localStorage.clear();
  useBridgeStore.getState().reset();
});

describe('LocalScope', () => {
  it('offers a connect card when the bridge is not connected', () => {
    renderWithProviders(
      <LocalScope>
        <Probe />
      </LocalScope>,
      { config: serverConfig },
    );
    expect(screen.getByText(/not connected to a local refgenie/i)).toBeInTheDocument();
    expect(screen.queryByTestId('base')).toBeNull();
  });

  it('points the clients at the bridge target and names the instance', () => {
    useBridgeStore
      .getState()
      .setConnected(permissivePing, 8080, 'http://localhost:8080');
    renderWithProviders(
      <LocalScope>
        <Probe />
      </LocalScope>,
      { config: serverConfig },
    );

    expect(screen.getByTestId('base')).toHaveTextContent('http://localhost:8080/v4');
    expect(screen.getByTestId('route-base')).toHaveTextContent('/local');
    expect(screen.getByTestId('service')).toHaveTextContent('my laptop refgenie');
    // The banner is the segregation: local data is never rendered unlabelled.
    const banner = screen.getByText(/lives on your computer/i);
    expect(banner).toHaveTextContent(/my laptop refgenie/);
    expect(banner.closest('.rg-banner')).not.toBeNull();
  });

  it('masks every command capability, whatever the ping claimed', () => {
    useBridgeStore
      .getState()
      .setConnected(permissivePing, 8080, 'http://localhost:8080');
    renderWithProviders(
      <LocalScope>
        <Probe />
      </LocalScope>,
      { config: serverConfig },
    );

    const on = screen.getByTestId('caps').textContent?.split(',').filter(Boolean) ?? [];
    expect(on.sort()).toEqual([...ALLOWED].sort());
    for (const key of CAPABILITY_KEYS) {
      if (!ALLOWED.includes(key)) expect(on).not.toContain(key);
    }
  });
});
