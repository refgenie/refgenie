import { beforeEach, describe, expect, it } from 'vitest';
import { act, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { BridgeJobsDock } from './BridgeJobsDock';
import { BridgeControl } from './BridgeControl';
import { server } from '../../test/server';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import { makeJob } from '../../test/jobFixtures';
import { useBridgeStore } from '../../stores/bridgeStore';
import { useBridgeJobStore } from '../../stores/bridgeJobStore';
import { coerceCapabilities } from '../../services/capabilities';
import type { PingResponse } from '../../services/bridge/contract';

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
  capabilities: coerceCapabilities({ pull: true, jobs: true }),
  bridge: { actions_cross_origin: true },
};

function track(id: string, instanceId = 'instance-1') {
  useBridgeJobStore.getState().track({
    id,
    instanceId,
    label: 'fasta/default',
    genomeDigest: 'g1',
    assetGroupName: 'fasta',
    assetName: 'default',
    submittedAt: '2026-10-27T00:00:00Z',
  });
}

/** Answers `statuses` in order, then repeats the last; counts every read. */
function serveJob(id: string, statuses: Array<'queued' | 'running' | 'succeeded' | 'failed'>) {
  const reads = { count: 0 };
  server.use(
    http.get(`${BRIDGE_ORIGIN}/v1/jobs/${id}`, () => {
      const status = statuses[Math.min(reads.count, statuses.length - 1)];
      reads.count += 1;
      return HttpResponse.json(
        makeJob({
          id,
          status,
          progress: status === 'running' ? { phase: 'download', percent: 40, message: null, bytes_done: null, bytes_total: null } : null,
          result:
            status === 'succeeded'
              ? { asset_digest: 'asset-digest-1', registry_path: 'hg38/fasta:default' }
              : null,
          error:
            status === 'failed'
              ? { code: 'pull_failed', message: 'disk full', detail: null, field: null }
              : null,
        }),
      );
    }),
  );
  return reads;
}

const sleep = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

beforeEach(() => {
  localStorage.clear();
  useBridgeStore.getState().reset();
  useBridgeJobStore.getState().clear();
});

describe('BridgeJobsDock', () => {
  it('renders nothing and polls nothing while disconnected', async () => {
    track('job-1');
    const reads = serveJob('job-1', ['running']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    await sleep(50);
    expect(screen.queryByRole('complementary', { name: /local pulls/i })).toBeNull();
    expect(reads.count).toBe(0);
  });

  it('reports a tracked pull to completion with a link into /local', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    track('job-1');
    serveJob('job-1', ['running', 'succeeded']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });

    expect(await screen.findByRole('progressbar')).toBeInTheDocument();
    const link = await screen.findByRole(
      'link',
      { name: /view hg38\/fasta:default/i },
      { timeout: 3000 },
    );
    expect(link).toHaveAttribute('href', '/local/assets/asset-digest-1');
    // The completion toast fires from the poller, not from a pull button.
    expect(screen.getByText(/fasta\/default is now on your local refgenie/i)).toBeInTheDocument();
  });

  it('shows the failure message from the job record', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    track('job-1');
    serveJob('job-1', ['failed']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    expect((await screen.findAllByText('disk full')).length).toBeGreaterThan(0);
    expect(screen.getByText('Failed')).toBeInTheDocument();
  });

  it('links a not-subscribed failure to the local prefilled subscribe form', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    useBridgeJobStore.getState().track({
      id: 'job-1',
      instanceId: 'instance-1',
      label: 'fasta/default',
      genomeDigest: 'g1',
      assetGroupName: 'fasta',
      assetName: 'default',
      serverUrl: 'https://example.org',
      submittedAt: '2026-10-27T00:00:00Z',
    });
    server.use(
      http.get(`${BRIDGE_ORIGIN}/v1/jobs/job-1`, () =>
        HttpResponse.json(
          makeJob({
            id: 'job-1',
            status: 'failed',
            error: {
              code: 'server_not_subscribed',
              message: 'not subscribed',
              detail: null,
              field: null,
            },
          }),
        ),
      ),
    );
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    const link = await screen.findByRole('link', { name: /subscribe on your local refgenie/i });
    expect(link).toHaveAttribute(
      'href',
      `${BRIDGE_ORIGIN}/servers?subscribe=${encodeURIComponent('https://example.org')}`,
    );
  });

  it('offers no subscribe link for other failures', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    track('job-1');
    serveJob('job-1', ['failed']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    expect((await screen.findAllByText('disk full')).length).toBeGreaterThan(0);
    expect(screen.queryByRole('link', { name: /subscribe on your local refgenie/i })).toBeNull();
  });

  it('marks a job the local refgenie no longer knows as lost and stops polling it', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    track('job-1');
    let reads = 0;
    server.use(
      http.get(`${BRIDGE_ORIGIN}/v1/jobs/job-1`, () => {
        reads += 1;
        return HttpResponse.json(
          { error: { code: 'not_found', message: 'no such job' } },
          { status: 404 },
        );
      }),
    );
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    expect(await screen.findByText(/no longer has a record/i)).toBeInTheDocument();
    await sleep(1200);
    expect(reads).toBe(1);
  });

  it('stops polling when the bridge disconnects', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    track('job-1');
    const reads = serveJob('job-1', ['running']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    await waitFor(() => expect(reads.count).toBe(1));

    act(() => useBridgeStore.getState().setAbsent(false));
    const after = reads.count;
    await sleep(1200);
    expect(reads.count).toBe(after);
    // The job is still remembered, so a reconnect resumes it.
    expect(useBridgeJobStore.getState().entries).toHaveLength(1);
  });

  it('ignores jobs tracked on a different local instance', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    track('other', 'instance-2');
    const reads = serveJob('other', ['running']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    await sleep(50);
    expect(screen.queryByRole('complementary', { name: /local pulls/i })).toBeNull();
    expect(reads.count).toBe(0);
  });

  it('forgets tracked jobs on an explicit Disconnect', async () => {
    useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
    track('job-1');
    serveJob('job-1', ['running']);
    renderWithProviders(
      <>
        <BridgeControl />
        <BridgeJobsDock />
      </>,
      { config: serverConfig },
    );
    expect(await screen.findByRole('complementary', { name: /local pulls/i })).toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: /disconnect/i }));
    expect(useBridgeJobStore.getState().entries).toEqual([]);
    expect(localStorage.getItem('refgenie.localBridge.jobs')).toBeNull();
    expect(screen.queryByRole('complementary', { name: /local pulls/i })).toBeNull();
  });
});
