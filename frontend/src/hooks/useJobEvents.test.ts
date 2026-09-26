/**
 * Transport tests for the job event stream. Every one of them runs without a
 * backend: the `EventSource` is a fake and the polling fallback goes through
 * the `fetchJobs` seam.
 */

import { act, renderHook } from '@testing-library/react';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { useJobEvents } from './useJobEvents';
import { useJobStore } from '../stores/jobStore';
import { createWrapper } from '../test/renderWithProviders';
import {
  FakeEventSource,
  createdSources,
  fakeEventSourceFactory,
  resetFakeEventSources,
} from '../test/FakeEventSource';
import { makeJob, page } from '../test/jobFixtures';
import type { JobsFetcher } from '../services/jobs';

const HEARTBEAT_MS = 1000;

function setup(options: Parameters<typeof useJobEvents>[0] = {}) {
  const fetchJobs = vi.fn(async () => page([])) as unknown as JobsFetcher;
  const rendered = renderHook(
    () =>
      useJobEvents({
        eventSourceFactory: fakeEventSourceFactory,
        fetchJobs,
        heartbeatMs: HEARTBEAT_MS,
        backoffMs: [1000, 2000, 4000],
        failoverAfter: 3,
        pollIntervalMs: { active: 2000, idle: 15000 },
        sseRetryMs: 60_000,
        ...options,
      }),
    { wrapper: createWrapper() },
  );
  return { fetchJobs: fetchJobs as unknown as ReturnType<typeof vi.fn>, rendered };
}

/** Flush the promise chain that `hydrate()` sits in front of `openSource()`. */
async function flush() {
  await act(async () => {
    await Promise.resolve();
    await Promise.resolve();
    await Promise.resolve();
    await Promise.resolve();
  });
}

async function advance(ms: number) {
  await act(async () => {
    await vi.advanceTimersByTimeAsync(ms);
  });
}

function latest(): FakeEventSource {
  return createdSources[createdSources.length - 1];
}

beforeEach(() => {
  vi.useFakeTimers();
  resetFakeEventSources();
  useJobStore.getState().reset();
});

afterEach(() => {
  vi.useRealTimers();
  vi.restoreAllMocks();
});

describe('useJobEvents', () => {
  it('hydrates and opens the multiplexed stream on mount', async () => {
    const { fetchJobs } = setup();
    await flush();

    expect(fetchJobs.mock.calls[0][0]).toMatchObject({ status: 'active' });
    expect(fetchJobs.mock.calls[1][0]).toMatchObject({ status: 'terminal' });
    expect(createdSources).toHaveLength(1);
    expect(createdSources[0].url).toContain('/v1/jobs/events');
  });

  it('applies progress, log and done events and invalidates on success', async () => {
    setup();
    await flush();
    const source = latest();
    act(() => source.open());

    act(() => {
      source.emit('status', { seq: 1, job_id: 'j-1', job: makeJob({ status: 'running' }) });
      source.emit('progress', {
        seq: 2,
        job_id: 'j-1',
        phase: 'download',
        percent: 38.2,
        message: 'hg38/fasta',
      });
      source.emit('log', { seq: 3, job_id: 'j-1', line: 'downloading…' });
    });

    let state = useJobStore.getState();
    expect(state.jobs['j-1'].progress?.percent).toBeCloseTo(38.2);
    expect(state.logs['j-1'].map((line) => line.text)).toEqual(['downloading…']);

    act(() => {
      source.emit('done', {
        seq: 4,
        job_id: 'j-1',
        job: makeJob({ status: 'succeeded', finished_at: '2026-08-12T10:05:00Z' }),
      });
    });

    state = useJobStore.getState();
    expect(state.jobs['j-1'].status).toBe('succeeded');
    expect(state.lastSeq).toBe(4);
  });

  it('ignores a replayed event whose seq is at or below the cursor', async () => {
    setup();
    await flush();
    const source = latest();
    act(() => source.open());

    const frame = { seq: 7, job_id: 'j-1', line: 'once' };
    act(() => {
      source.emit('log', frame);
      source.emit('log', frame);
    });

    expect(useJobStore.getState().logs['j-1']).toHaveLength(1);

    act(() => source.emit('log', { seq: 8, job_id: 'j-1', line: 'twice' }));
    expect(useJobStore.getState().logs['j-1']).toHaveLength(2);
  });

  it('reconnects after one failure, resuming from the cursor', async () => {
    setup();
    await flush();
    const source = latest();
    act(() => source.open());
    act(() => source.emit('heartbeat', { seq: 5 }));

    act(() => source.fail());
    expect(source.closed).toBe(true);
    expect(createdSources).toHaveLength(1);

    await advance(1000);
    expect(createdSources).toHaveLength(2);
    expect(latest().url).toContain('since=5');
  });

  it('fails over to polling after three consecutive failures', async () => {
    const { fetchJobs } = setup();
    await flush();

    act(() => latest().fail());
    await advance(1000);
    act(() => latest().fail());
    await advance(2000);
    act(() => latest().fail());
    await flush();

    expect(createdSources).toHaveLength(3);
    expect(useJobStore.getState().connection).toBe('polling');

    const before = fetchJobs.mock.calls.length;
    await advance(15_000);
    expect(fetchJobs.mock.calls.length).toBeGreaterThan(before);
    // Still three: polling means we stopped reopening.
    expect(createdSources).toHaveLength(3);
  });

  it('returns to the live stream when the SSE retry succeeds', async () => {
    // A long heartbeat window so the silence watchdog does not tear the
    // recovered connection down again inside the assertion window.
    const { fetchJobs } = setup({ heartbeatMs: 100_000 });
    await flush();
    act(() => latest().fail());
    await advance(1000);
    act(() => latest().fail());
    await advance(2000);
    act(() => latest().fail());
    await flush();
    expect(useJobStore.getState().connection).toBe('polling');

    await advance(60_000);
    expect(createdSources).toHaveLength(4);
    act(() => latest().open());

    expect(useJobStore.getState().connection).toBe('open');
    const after = fetchJobs.mock.calls.length;
    await advance(60_000);
    expect(fetchJobs.mock.calls.length).toBe(after);
  });

  it('pauses in a hidden tab and re-hydrates when it comes back', async () => {
    const { fetchJobs } = setup();
    await flush();
    act(() => latest().open());
    const opened = createdSources.length;

    Object.defineProperty(document, 'visibilityState', {
      configurable: true,
      get: () => 'hidden',
    });
    await act(async () => {
      document.dispatchEvent(new Event('visibilitychange'));
    });
    expect(createdSources[opened - 1].closed).toBe(true);
    expect(useJobStore.getState().connection).toBe('offline');

    const before = fetchJobs.mock.calls.length;
    Object.defineProperty(document, 'visibilityState', {
      configurable: true,
      get: () => 'visible',
    });
    await act(async () => {
      document.dispatchEvent(new Event('visibilitychange'));
    });
    await flush();

    expect(fetchJobs.mock.calls.length).toBeGreaterThan(before);
    expect(createdSources.length).toBe(opened + 1);
  });

  it('skips a malformed frame and keeps the stream alive', async () => {
    vi.spyOn(console, 'warn').mockImplementation(() => {});
    setup();
    await flush();
    const source = latest();
    act(() => source.open());

    act(() => source.emitRaw('progress', 'not json{'));
    expect(source.closed).toBe(false);

    act(() => source.emit('log', { seq: 2, job_id: 'j-1', line: 'still here' }));
    expect(useJobStore.getState().logs['j-1']).toHaveLength(1);
  });

  it('reconnects when three heartbeat intervals pass in silence', async () => {
    setup();
    await flush();
    act(() => latest().open());
    expect(createdSources).toHaveLength(1);

    // Silence, not an error: the browser would never fire onerror here.
    await advance(HEARTBEAT_MS * 3);
    await advance(1000);
    expect(createdSources.length).toBeGreaterThan(1);
  });

  it('closes the source and leaves no timers on unmount', async () => {
    const { rendered } = setup();
    await flush();
    act(() => latest().open());

    act(() => rendered.unmount());

    expect(latest().closed).toBe(true);
    expect(vi.getTimerCount()).toBe(0);
  });
});
