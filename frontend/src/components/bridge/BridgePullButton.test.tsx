import { beforeEach, describe, expect, it } from 'vitest';
import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { BridgePullButton } from './BridgePullButton';
import { server } from '../../test/server';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import { useBridgeStore } from '../../stores/bridgeStore';
import { coerceCapabilities } from '../../services/capabilities';
import type { PingResponse } from '../../services/bridge/contract';
import type { UiConfig } from '../../types/ui';

const BRIDGE_ORIGIN = 'http://localhost:8080';
const GENOME = '0000000000000000000000000000cafe';

/**
 * The page's own backend is cross-origin here, which is the real deployment:
 * refgenie.org's SPA reading api.refgenie.org.
 */
const remoteConfig: UiConfig = {
  ...serverConfig,
  api_base: 'https://api.refgenie.org/v4',
};

function ping(overrides: Partial<PingResponse> = {}): PingResponse {
  return {
    service: 'refgenie',
    bridge_version: 1,
    mode: 'local',
    refgenie_version: '1.0.0',
    instance_id: 'instance-1',
    instance_label: 'local refgenie',
    bridge_mode: 'full',
    action_header: 'X-Refgenie-Action',
    capabilities: coerceCapabilities({ pull: true }),
    bridge: { actions_cross_origin: true },
    ...overrides,
  };
}

function connect(p: PingResponse) {
  useBridgeStore.getState().setConnected(p, 8080, BRIDGE_ORIGIN);
}

function jobRecord(status: string) {
  return {
    id: 'job-1',
    kind: 'pull',
    status,
    queue_position: null,
    label: 'pull',
    target: {
      genome_digest: GENOME,
      genome_name: null,
      asset_group_name: 'fasta',
      asset_name: 'default',
    },
    progress: status === 'running' ? { phase: 'download', percent: 40 } : null,
    created_at: '2026-08-26T00:00:00Z',
    started_at: null,
    finished_at: null,
    result: null,
    error: null,
    log_lines: 0,
  };
}

beforeEach(() => {
  localStorage.clear();
  useBridgeStore.getState().reset();
});

describe('BridgePullButton', () => {
  it('gates on the PING capability, not the page capability', () => {
    // The page's own `pull` is true, the local instance's is false: nothing.
    connect(ping({ capabilities: coerceCapabilities({ pull: false }) }));
    const { unmount } = renderWithProviders(
      <BridgePullButton genomeDigest={GENOME} assetGroupName="fasta" />,
      { config: { ...remoteConfig, capabilities: coerceCapabilities({ pull: true }) } },
    );
    expect(screen.queryByRole('button', { name: /pull to my refgenie/i })).toBeNull();
    unmount();

    // And the reverse: the page cannot pull (it is refgenie.org), the local
    // instance can. This is the regression that would kill the feature.
    connect(ping());
    renderWithProviders(
      <BridgePullButton genomeDigest={GENOME} assetGroupName="fasta" />,
      { config: remoteConfig },
    );
    expect(
      screen.getByRole('button', { name: /pull to my refgenie/i }),
    ).toBeInTheDocument();
  });

  it('renders nothing when the bridge is not connected', () => {
    renderWithProviders(
      <BridgePullButton genomeDigest={GENOME} assetGroupName="fasta" />,
      { config: remoteConfig },
    );
    expect(screen.queryByRole('button', { name: /pull to my refgenie/i })).toBeNull();
  });

  it('offers the deep link instead when cross-origin actions are off', () => {
    connect(ping({ bridge_mode: 'read', bridge: { actions_cross_origin: false } }));
    renderWithProviders(
      <BridgePullButton
        genomeDigest={GENOME}
        assetGroupName="fasta"
        assetName="default"
      />,
      { config: remoteConfig },
    );

    const button = screen.getByRole('button', { name: /pull to my refgenie/i });
    expect(button).toBeDisabled();
    expect(button).toHaveAttribute(
      'title',
      expect.stringContaining('refgenie dash --bridge full'),
    );

    const link = screen.getByRole('link', { name: /open in local refgenie/i });
    expect(link).toHaveAttribute(
      'href',
      `${BRIDGE_ORIGIN}/pull?server=${encodeURIComponent('https://api.refgenie.org')}` +
        `&genome=${GENOME}&asset_group=fasta&asset=default`,
    );
  });

  it('submits force:false with the action header, then polls to success', async () => {
    connect(ping());

    let body: Record<string, unknown> | undefined;
    let actionHeader: string | null = null;
    let jobReads = 0;
    server.use(
      http.post(`${BRIDGE_ORIGIN}/v1/actions/pull`, async ({ request }) => {
        body = (await request.json()) as Record<string, unknown>;
        actionHeader = request.headers.get('x-refgenie-action');
        return HttpResponse.json(
          { job_id: 'job-1', kind: 'pull', status: 'queued', duplicate: false },
          { status: 202 },
        );
      }),
      http.get(`${BRIDGE_ORIGIN}/v1/jobs/:id`, () => {
        jobReads += 1;
        return HttpResponse.json(jobRecord(jobReads === 1 ? 'running' : 'succeeded'));
      }),
    );

    renderWithProviders(
      <BridgePullButton
        genomeDigest={GENOME}
        assetGroupName="fasta"
        assetName="default"
      />,
      { config: remoteConfig },
    );

    await userEvent.click(screen.getByRole('button', { name: /pull to my refgenie/i }));

    await waitFor(() => expect(body).toBeDefined());
    expect(actionHeader).toBe('1');
    expect(body).toEqual({
      asset_group: 'fasta',
      genome_digest: GENOME,
      asset: 'default',
      server_url: 'https://api.refgenie.org',
      force: false,
    });
    // Exactly one genome reference: `genome` is omitted, not nulled.
    expect(body).not.toHaveProperty('genome');

    // The first poll is not terminal, so the inline progress bar renders.
    expect(await screen.findByRole('progressbar')).toBeInTheDocument();

    await waitFor(
      () => expect(screen.getByText(/now on your local refgenie/i)).toBeInTheDocument(),
      { timeout: 4000 },
    );
  });

  it('falls back to the deep link when the local side is not subscribed', async () => {
    connect(ping());
    const message =
      'This refgenie is not subscribed to that server, so a remote page cannot pull from it.';
    server.use(
      http.post(`${BRIDGE_ORIGIN}/v1/actions/pull`, () =>
        HttpResponse.json(
          { error: { code: 'server_not_subscribed', message } },
          { status: 403 },
        ),
      ),
    );

    renderWithProviders(
      <BridgePullButton
        genomeDigest={GENOME}
        assetGroupName="fasta"
        assetName="default"
      />,
      { config: remoteConfig },
    );
    // Cross-origin actions are on, so there is no fallback link yet.
    expect(screen.queryByRole('link', { name: /open in local refgenie/i })).toBeNull();

    await userEvent.click(screen.getByRole('button', { name: /pull to my refgenie/i }));

    expect(await screen.findByText(message)).toBeInTheDocument();
    const link = await screen.findByRole('link', { name: /open in local refgenie/i });
    expect(link).toHaveAttribute(
      'href',
      `${BRIDGE_ORIGIN}/pull?server=${encodeURIComponent('https://api.refgenie.org')}` +
        `&genome=${GENOME}&asset_group=fasta&asset=default`,
    );
  });
});
