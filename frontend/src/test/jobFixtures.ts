import type { Job, JobKind, JobStatus } from '../services/contracts';
import type { Paginated } from '../types/pagination';

export function makeJob(overrides: Partial<Job> = {}): Job {
  const kind: JobKind = overrides.kind ?? 'pull';
  const status: JobStatus = overrides.status ?? 'running';
  return {
    id: 'j-1',
    kind,
    status,
    queue_position: null,
    label: 'pull hg38/fasta:default from refgenomes.databio.org',
    target: {
      genome_digest: 'genome-digest-1',
      genome_name: 'hg38',
      asset_group_name: 'fasta',
      asset_name: 'default',
    },
    progress: null,
    created_at: '2026-08-12T10:00:00Z',
    started_at: status === 'queued' ? null : '2026-08-12T10:00:01Z',
    finished_at: null,
    result: null,
    error: null,
    log_lines: 0,
    cancellable: true,
    ...overrides,
  };
}

export function page<T>(items: T[]): Paginated<T> {
  return { items, pagination: { offset: 0, limit: 50, total: items.length } };
}
