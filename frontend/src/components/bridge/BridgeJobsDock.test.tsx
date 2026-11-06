import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { BridgeJobsDock } from './BridgeJobsDock';
import { BridgeControl } from './BridgeControl';
import { server } from '../../test/server';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import { BRIDGE_ORIGIN, makePing } from '../../test/bridgeFixtures';
import { makeJob } from '../../test/jobFixtures';
import { useBridgeStore } from '../../stores/bridgeStore';
import { useBridgeJobStore } from '../../stores/bridgeJobStore';
import { coerceCapabilities } from '../../services/capabilities';
import { BRIDGE_POLL_INTERVAL_MS } from '../../services/bridge/pull';
import type { JobError } from '../../services/contracts';

const ping = makePing({ capabilities: coerceCapabilities({ pull: true, jobs: true }) });

const DISK_FULL: JobError = { code: 'pull_failed', message: 'disk full', detail: null, field: null };

function connect() {
  useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
}

function track(id: string, instanceId = 'instance-1') {
  useBridgeJobStore.getState().track({
    id,
    instanceId,
    label: 'fasta/default',
    genomeDigest: 'g1',
    assetGroupName: 'fasta',
    assetName: 'default',
    serverUrl: 'https://example.org',
    submittedAt: '2026-10-27T00:00:00Z',
  });
}

/** Answers `statuses` in order, then repeats the last; counts every read. */
function serveJob(
  id: string,
  statuses: Array<'queued' | 'running' | 'succeeded' | 'failed'>,
  error: JobError = DISK_FULL,
) {
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
          error: status === 'failed' ? error : null,
        }),
      );
    }),
  );
  return reads;
}

/** Runs the poller past `polls` intervals in virtual time. */
async function advancePolls(polls = 1) {
  await act(() => vi.advanceTimersByTimeAsync(polls * BRIDGE_POLL_INTERVAL_MS + 1));
}

beforeEach(() => {
  // Real time still flows so findBy* and MSW work; polls are driven by advancePolls.
  vi.useFakeTimers({ shouldAdvanceTime: true });
  localStorage.clear();
  useBridgeStore.getState().reset();
  useBridgeJobStore.getState().clear();
});

afterEach(() => {
  vi.useRealTimers();
});

describe('BridgeJobsDock', () => {
  it('renders nothing and polls nothing while disconnected', async () => {
    track('job-1');
    const reads = serveJob('job-1', ['running']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    await advancePolls();
    expect(screen.queryByRole('complementary', { name: /local pulls/i })).toBeNull();
    expect(reads.count).toBe(0);
  });

  it('reports a tracked pull to completion with a link into /local', async () => {
    connect();
    track('job-1');
    serveJob('job-1', ['running', 'succeeded']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });

    expect(await screen.findByRole('progressbar')).toBeInTheDocument();
    await advancePolls();
    const link = await screen.findByRole('link', { name: /view hg38\/fasta:default/i });
    expect(link).toHaveAttribute('href', '/local/assets/asset-digest-1');
    // The completion toast fires from the poller, not from a pull button.
    expect(screen.getByText(/fasta\/default is now on your local refgenie/i)).toBeInTheDocument();
  });

  it.each([
    { name: 'an ordinary failure', error: DISK_FULL, subscribeHref: null },
    {
      name: 'a not-subscribed failure',
      error: { code: 'server_not_subscribed', message: 'not subscribed', detail: null, field: null },
      subscribeHref: `${BRIDGE_ORIGIN}/servers?subscribe=${encodeURIComponent('https://example.org')}`,
    },
  ])('shows $name with the right subscribe link', async ({ error, subscribeHref }) => {
    connect();
    track('job-1');
    serveJob('job-1', ['failed'], error);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });

    expect((await screen.findAllByText(error.message)).length).toBeGreaterThan(0);
    expect(screen.getByText('Failed')).toBeInTheDocument();
    const link = screen.queryByRole('link', { name: /subscribe on your local refgenie/i });
    if (subscribeHref) {
      expect(link).toHaveAttribute('href', subscribeHref);
    } else {
      expect(link).toBeNull();
    }
  });

  it('marks a job the local refgenie no longer knows as lost and stops polling it', async () => {
    connect();
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
    await advancePolls(2);
    expect(reads).toBe(1);
  });

  it('stops polling when the bridge disconnects', async () => {
    connect();
    track('job-1');
    const reads = serveJob('job-1', ['running']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    expect(await screen.findByRole('progressbar')).toBeInTheDocument();
    expect(reads.count).toBe(1);

    act(() => useBridgeStore.getState().setAbsent(false));
    await advancePolls(2);
    expect(reads.count).toBe(1);
    // The job is still remembered, so a reconnect resumes it.
    expect(useBridgeJobStore.getState().entries).toHaveLength(1);
  });

  it('ignores jobs tracked on a different local instance', async () => {
    connect();
    track('other', 'instance-2');
    const reads = serveJob('other', ['running']);
    renderWithProviders(<BridgeJobsDock />, { config: serverConfig });
    await advancePolls();
    expect(screen.queryByRole('complementary', { name: /local pulls/i })).toBeNull();
    expect(reads.count).toBe(0);
  });

  it('forgets tracked jobs on an explicit Disconnect', async () => {
    connect();
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
