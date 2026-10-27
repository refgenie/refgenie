import { beforeEach, describe, expect, it } from 'vitest';
import {
  LOG_RING_CAP,
  selectQueued,
  targetKey,
  useJobStore,
} from './jobStore';
import { makeJob } from '../test/jobFixtures';

beforeEach(() => useJobStore.getState().reset());

describe('targetKey', () => {
  it('prefers the digest, because an alias and a digest are the same work', () => {
    expect(
      targetKey('pull', {
        genome_digest: 'abc',
        genome_name: 'hg38',
        asset_group_name: 'fasta',
        asset_name: 'default',
      }),
    ).toBe('pull:abc/fasta:default');
  });

  it('falls back to the name, and stands in a wildcard for a missing asset', () => {
    expect(
      targetKey('build', {
        genome_digest: null,
        genome_name: 'hg38',
        asset_group_name: 'bowtie2_index',
        asset_name: null,
      }),
    ).toBe('build:hg38/bowtie2_index:*');
  });
});

describe('activeJobFor', () => {
  it('matches a running job and a queued one, and stops at terminal', () => {
    const store = useJobStore.getState();
    const job = makeJob({ id: 'j-1', status: 'running' });
    store.upsertJob(job);
    const key = targetKey(job.kind, job.target);

    expect(useJobStore.getState().activeJobFor(key)?.id).toBe('j-1');

    useJobStore.getState().upsertJob({ ...job, status: 'queued' });
    expect(useJobStore.getState().activeJobFor(key)?.id).toBe('j-1');

    useJobStore.getState().upsertJob({ ...job, status: 'succeeded' });
    expect(useJobStore.getState().activeJobFor(key)).toBeNull();
  });
});

describe('selectQueued', () => {
  it('sorts by queue position within a kind, since the queues are independent', () => {
    const jobs = {
      a: makeJob({ id: 'a', kind: 'pull', status: 'queued', queue_position: 1 }),
      b: makeJob({ id: 'b', kind: 'pull', status: 'queued', queue_position: 0 }),
      c: makeJob({ id: 'c', kind: 'build', status: 'queued', queue_position: 0 }),
    };
    expect(selectQueued(jobs).map((job) => job.id)).toEqual(['c', 'b', 'a']);
  });
});

describe('log ring buffer', () => {
  it('caps at 500 lines and drops from the front', () => {
    const store = useJobStore.getState();
    const lines = Array.from({ length: LOG_RING_CAP + 10 }, (_, i) => `line ${i}`);
    store.appendLog('j-1', lines);

    const log = useJobStore.getState().logs['j-1'];
    expect(log).toHaveLength(LOG_RING_CAP);
    expect(log[0].text).toBe('line 10');
    expect(log[log.length - 1].text).toBe(`line ${LOG_RING_CAP + 9}`);
  });

  it('tags pipeline output so tool noise reads differently from phase messages', () => {
    useJobStore.getState().appendLog('j-1', ['bowtie2-build ...'], 'pipeline');
    expect(useJobStore.getState().logs['j-1'][0].source).toBe('pipeline');
  });
});
