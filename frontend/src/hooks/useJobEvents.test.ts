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

function setup() {
  const fetchJobs = vi.fn(async () => page([])) as unknown as JobsFetcher;
  const rendered = renderHook(
    () => useJobEvents({ eventSourceFactory: fakeEventSourceFactory, fetchJobs }),
    { wrapper: createWrapper() },
  );
  return { fetchJobs: fetchJobs as unknown as ReturnType<typeof vi.fn>, rendered };
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
    await advance(0);

    expect(fetchJobs.mock.calls[0][0]).toMatchObject({ status: 'active' });
    expect(fetchJobs.mock.calls[1][0]).toMatchObject({ status: 'terminal' });
    expect(createdSources).toHaveLength(1);
    expect(createdSources[0].url).toContain('/v1/jobs/events');
  });

  it('applies events, skipping malformed frames', async () => {
    vi.spyOn(console, 'warn').mockImplementation(() => {});
    setup();
    await advance(0);
    const source = latest();
    act(() => source.open());
    expect(useJobStore.getState().connection).toBe('open');

    act(() => {
      source.emit('status', {
        seq: 1,
        job_id: 'j-1',
        job: makeJob({ status: 'running' }),
      });
      source.emitRaw('progress', 'not json{');
      source.emit('progress', {
        seq: 2,
        job_id: 'j-1',
        phase: 'download',
        percent: 38.2,
      });
      source.emit('log', { seq: 3, job_id: 'j-1', line: 'downloading…' });
    });
    let state = useJobStore.getState();
    expect(state.jobs['j-1'].progress?.percent).toBeCloseTo(38.2);
    expect(state.logs['j-1'].map((line) => line.text)).toEqual(['downloading…']);

    act(() =>
      source.emit('done', {
        seq: 4,
        job_id: 'j-1',
        job: makeJob({ status: 'succeeded' }),
      }),
    );
    state = useJobStore.getState();
    expect(state.jobs['j-1'].status).toBe('succeeded');
    expect(state.lastSeq).toBe(4);
  });

  it('on failure, re-reads the job list and reopens from the cursor', async () => {
    const { fetchJobs } = setup();
    await advance(0);
    const source = latest();
    act(() => source.open());
    act(() => source.emit('heartbeat', { seq: 5 }));

    act(() => source.fail());
    expect(source.closed).toBe(true);
    expect(useJobStore.getState().connection).toBe('polling');

    const before = fetchJobs.mock.calls.length;
    await advance(2000);
    expect(fetchJobs.mock.calls.length).toBeGreaterThan(before);
    expect(createdSources).toHaveLength(2);
    expect(latest().url).toContain('since=5');
  });

  it('closes the source and leaves no timers on unmount', async () => {
    const { rendered } = setup();
    await advance(0);
    act(() => latest().fail());

    act(() => rendered.unmount());

    expect(latest().closed).toBe(true);
    expect(vi.getTimerCount()).toBe(0);
  });
});
