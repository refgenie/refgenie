import { beforeEach, describe, expect, it, vi } from 'vitest';
import { BRIDGE_JOB_CAP, isActiveEntry, useBridgeJobStore } from './bridgeJobStore';
import { makeJob } from '../test/jobFixtures';

const KEY = 'refgenie.localBridge.jobs';

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

beforeEach(() => {
  localStorage.clear();
  useBridgeJobStore.getState().clear();
});

describe('bridgeJobStore', () => {
  it('writes every change through to localStorage', () => {
    track('job-1');
    useBridgeJobStore.getState().update('job-1', makeJob({ id: 'job-1', status: 'running' }));
    const stored = JSON.parse(localStorage.getItem(KEY) ?? '[]');
    expect(stored).toHaveLength(1);
    expect(stored[0]).toMatchObject({ id: 'job-1', instanceId: 'instance-1' });
    expect(stored[0].job.status).toBe('running');

    useBridgeJobStore.getState().dismiss('job-1');
    expect(localStorage.getItem(KEY)).toBeNull();
  });

  it('resumes remembered jobs on a fresh page load and drops malformed ones', async () => {
    localStorage.setItem(
      KEY,
      JSON.stringify([
        {
          id: 'job-9',
          instanceId: 'instance-1',
          label: 'bowtie2_index',
          genomeDigest: 'g1',
          assetGroupName: 'bowtie2_index',
          assetName: null,
          submittedAt: '2026-10-27T00:00:00Z',
          job: null,
          lost: false,
        },
        { id: 7 },
        'junk',
      ]),
    );
    vi.resetModules();
    const fresh = await import('./bridgeJobStore');
    const entries = fresh.useBridgeJobStore.getState().entries;
    expect(entries.map((entry) => entry.id)).toEqual(['job-9']);
    expect(fresh.isActiveEntry(entries[0])).toBe(true);
  });

  it('reads unparseable storage as no jobs', async () => {
    localStorage.setItem(KEY, '{not json');
    vi.resetModules();
    const fresh = await import('./bridgeJobStore');
    expect(fresh.useBridgeJobStore.getState().entries).toEqual([]);
  });

  it('keeps one card per job id, newest first, capped', () => {
    for (let i = 0; i < BRIDGE_JOB_CAP + 5; i += 1) track(`job-${i}`);
    track('job-3');
    const ids = useBridgeJobStore.getState().entries.map((entry) => entry.id);
    expect(ids).toHaveLength(BRIDGE_JOB_CAP);
    expect(ids[0]).toBe('job-3');
    expect(new Set(ids).size).toBe(ids.length);
  });

  it('clears finished and lost jobs of one instance only', () => {
    track('done');
    track('lost');
    track('running');
    track('elsewhere', 'instance-2');
    const store = useBridgeJobStore.getState();
    store.update('done', makeJob({ id: 'done', status: 'succeeded' }));
    store.update('running', makeJob({ id: 'running', status: 'running' }));
    store.update('elsewhere', makeJob({ id: 'elsewhere', status: 'failed' }));
    store.markLost('lost');

    const entries = useBridgeJobStore.getState().entries;
    expect(entries.filter(isActiveEntry).map((entry) => entry.id)).toEqual(['running']);

    useBridgeJobStore.getState().clearFinished('instance-1');
    expect(useBridgeJobStore.getState().entries.map((entry) => entry.id).sort()).toEqual([
      'elsewhere',
      'running',
    ]);
  });
});
