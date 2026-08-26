/**
 * Live job state.
 *
 * zustand rather than React context specifically because log lines arrive at
 * tens of events per second during a build: a context holding that state
 * re-renders every consumer, while a zustand selector lets the log tail
 * re-render alone.
 *
 * The store is a plain module with no React import beyond the hook factory, so
 * `useJobEvents` can drive it and tests can exercise it without rendering.
 */

import { create } from 'zustand';
import { INVALIDATION_FOR_JOB_KIND, invalidate } from '../services/invalidation';
import { isTerminal } from '../services/contracts';
import type {
  Job,
  JobEvent,
  JobKind,
  JobProgressInfo,
  JobStatus,
  JobTarget,
} from '../services/contracts';

/** Per-job log ring buffer cap. The full log comes from `/v1/jobs/{id}/log`. */
export const LOG_RING_CAP = 500;

/** How long a finished job stays on the console so a fast pull can be read. */
export const TERMINAL_LINGER_SECONDS = 60;

export type ConnectionState = 'connecting' | 'open' | 'polling' | 'offline';

export interface LogLine {
  text: string;
  /** `pipeline` is tee'd subprocess output; it renders visually distinct. */
  source: 'refgenie' | 'pipeline';
}

export interface JobStoreState {
  jobs: Record<string, Job>;
  logs: Record<string, LogLine[]>;
  /** Highest manager-global seq applied. Replay below this is ignored. */
  lastSeq: number;
  connection: ConnectionState;
  /** Job whose card the UI should scroll to and highlight, if any. */
  focusedJobId: string | null;

  upsertJob: (job: Job) => void;
  upsertMany: (jobs: readonly Job[]) => void;
  registerQueued: (job: Job) => void;
  applyEvent: (event: JobEvent) => void;
  appendLog: (jobId: string, lines: readonly string[], source?: LogLine['source']) => void;
  setConnection: (connection: ConnectionState) => void;
  setLastSeq: (seq: number) => void;
  focusJob: (jobId: string | null) => void;
  activeJobFor: (targetKey: string) => Job | null;
  clearTerminal: () => void;
  reset: () => void;
}

/**
 * The duplicate/active-job identity. `genome_digest` first because a pull
 * submitted by alias and one submitted by digest are the same work.
 */
export function targetKey(kind: JobKind, target: JobTarget | null | undefined): string {
  const genome = target?.genome_digest ?? target?.genome_name ?? '?';
  const group = target?.asset_group_name ?? '?';
  const name = target?.asset_name ?? '*';
  return `${kind}:${genome}/${group}:${name}`;
}

export function jobTargetKey(job: Job): string {
  return targetKey(job.kind, job.target);
}

function pushLines(
  existing: LogLine[] | undefined,
  lines: readonly string[],
  source: LogLine['source'],
): LogLine[] {
  const next = [...(existing ?? []), ...lines.map((text) => ({ text, source }))];
  return next.length > LOG_RING_CAP ? next.slice(next.length - LOG_RING_CAP) : next;
}

/** A record good enough to render before the first SSE frame arrives. */
export function provisionalJob(args: {
  id: string;
  kind: JobKind;
  label: string;
  target: JobTarget;
  status?: JobStatus;
  created_at?: string;
}): Job {
  return {
    id: args.id,
    kind: args.kind,
    status: args.status ?? 'queued',
    queue_position: null,
    label: args.label,
    target: args.target,
    progress: null,
    created_at: args.created_at ?? new Date().toISOString(),
    started_at: null,
    finished_at: null,
    result: null,
    error: null,
    log_lines: 0,
    cancellable: true,
  };
}

/** True when the payload looks like a whole job record rather than a patch. */
function looksLikeJob(value: unknown): value is Job {
  if (!value || typeof value !== 'object') return false;
  const candidate = value as Partial<Job>;
  return typeof candidate.id === 'string' && typeof candidate.status === 'string';
}

function mergeProgress(
  existing: JobProgressInfo | null,
  patch: Partial<JobProgressInfo>,
): JobProgressInfo {
  return {
    phase: patch.phase ?? existing?.phase ?? '',
    percent: patch.percent === undefined ? (existing?.percent ?? null) : patch.percent,
    message: patch.message === undefined ? (existing?.message ?? null) : patch.message,
    bytes_done:
      patch.bytes_done === undefined ? (existing?.bytes_done ?? null) : patch.bytes_done,
    bytes_total:
      patch.bytes_total === undefined ? (existing?.bytes_total ?? null) : patch.bytes_total,
  };
}

const EMPTY: Pick<JobStoreState, 'jobs' | 'logs' | 'lastSeq' | 'connection' | 'focusedJobId'> = {
  jobs: {},
  logs: {},
  lastSeq: 0,
  connection: 'connecting',
  focusedJobId: null,
};

export const useJobStore = create<JobStoreState>()((set, get) => ({
  ...EMPTY,

  upsertJob: (job) =>
    set((state) => ({ jobs: { ...state.jobs, [job.id]: { ...state.jobs[job.id], ...job } } })),

  upsertMany: (jobs) =>
    set((state) => {
      if (jobs.length === 0) return {};
      const next = { ...state.jobs };
      for (const job of jobs) next[job.id] = { ...next[job.id], ...job };
      return { jobs: next };
    }),

  // Not optimistic UI: the id came from the server, so the record is real.
  // Only the label is provisional, and the first `status` event replaces it.
  registerQueued: (job) =>
    set((state) => ({
      jobs: { ...state.jobs, [job.id]: { ...state.jobs[job.id], ...job } },
      focusedJobId: job.id,
    })),

  applyEvent: (event) => {
    const seq = typeof event.seq === 'number' ? event.seq : 0;
    // Idempotent replay: a reconnect resends from the cursor, and the server
    // is allowed to overlap. Anything at or below the cursor is already ours.
    if (seq > 0 && seq <= get().lastSeq) return;

    const jobId = event.job_id;

    if (event.type === 'heartbeat') {
      set({ lastSeq: Math.max(get().lastSeq, seq) });
      return;
    }

    if (!jobId) {
      set({ lastSeq: Math.max(get().lastSeq, seq) });
      return;
    }

    set((state) => {
      const existing = state.jobs[jobId];
      let jobs = state.jobs;
      let logs = state.logs;

      if (event.type === 'status') {
        const patch: Partial<Job> = looksLikeJob(event.job)
          ? event.job
          : {
              status: event.status,
              queue_position:
                event.queue_position === undefined ? null : event.queue_position,
            };
        jobs = { ...jobs, [jobId]: { ...existing, ...patch, id: jobId } as Job };
      } else if (event.type === 'progress') {
        const { seq: _seq, job_id: _jobId, type: _type, ...patch } = event;
        jobs = {
          ...jobs,
          [jobId]: {
            ...existing,
            id: jobId,
            progress: mergeProgress(existing?.progress ?? null, patch),
          } as Job,
        };
      } else if (event.type === 'log') {
        const lines = event.lines ?? (event.line === undefined ? [] : [event.line]);
        if (lines.length > 0) {
          const source = event.source === 'pipeline' ? 'pipeline' : 'refgenie';
          logs = { ...logs, [jobId]: pushLines(logs[jobId], lines, source) };
        }
      } else if (event.type === 'truncated') {
        const dropped = event.dropped ?? 0;
        logs = {
          ...logs,
          [jobId]: pushLines(
            logs[jobId],
            [`… ${dropped || 'some'} earlier lines dropped; open details for the full log`],
            'refgenie',
          ),
        };
      } else if (event.type === 'done') {
        const record = looksLikeJob(event.job)
          ? event.job
          : looksLikeJob(event)
            ? (event as unknown as Job)
            : undefined;
        if (record) {
          jobs = { ...jobs, [jobId]: { ...existing, ...record, id: jobId } };
        }
      }

      return { jobs, logs, lastSeq: Math.max(state.lastSeq, seq) };
    });

    if (event.type === 'done') {
      const finished = get().jobs[jobId];
      if (finished && finished.status === 'succeeded') {
        invalidate(INVALIDATION_FOR_JOB_KIND[finished.kind] ?? []);
      }
    }
  },

  appendLog: (jobId, lines, source = 'refgenie') =>
    set((state) => ({ logs: { ...state.logs, [jobId]: pushLines(state.logs[jobId], lines, source) } })),

  setConnection: (connection) => set({ connection }),

  setLastSeq: (seq) => set((state) => ({ lastSeq: Math.max(state.lastSeq, seq) })),

  focusJob: (focusedJobId) => set({ focusedJobId }),

  activeJobFor: (key) => {
    const match = Object.values(get().jobs).find(
      (job) => !isTerminal(job.status) && jobTargetKey(job) === key,
    );
    return match ?? null;
  },

  clearTerminal: () =>
    set((state) => {
      const jobs: Record<string, Job> = {};
      const logs: Record<string, LogLine[]> = {};
      for (const [id, job] of Object.entries(state.jobs)) {
        if (isTerminal(job.status)) continue;
        jobs[id] = job;
        if (state.logs[id]) logs[id] = state.logs[id];
      }
      return { jobs, logs, focusedJobId: null };
    }),

  reset: () => set({ ...EMPTY }),
}));

// === Pure selectors over the jobs map ===

function byCreatedAt(a: Job, b: Job): number {
  return a.created_at.localeCompare(b.created_at);
}

export function selectRunning(jobs: Record<string, Job>): Job[] {
  return Object.values(jobs).filter((job) => job.status === 'running').sort(byCreatedAt);
}

/**
 * Queued jobs, sorted by position WITHIN THEIR OWN KIND: pulls (2 slots) and
 * builds (1 slot) are independent queues, so a position is only meaningful
 * against its own executor.
 */
export function selectQueued(jobs: Record<string, Job>): Job[] {
  return Object.values(jobs)
    .filter((job) => job.status === 'queued')
    .sort((a, b) => {
      if (a.kind !== b.kind) return a.kind.localeCompare(b.kind);
      const posA = a.queue_position ?? Number.MAX_SAFE_INTEGER;
      const posB = b.queue_position ?? Number.MAX_SAFE_INTEGER;
      return posA !== posB ? posA - posB : byCreatedAt(a, b);
    });
}

/** Terminal jobs that finished within `seconds`, newest first. */
export function selectRecentTerminal(
  jobs: Record<string, Job>,
  seconds: number = TERMINAL_LINGER_SECONDS,
  now: number = Date.now(),
): Job[] {
  const cutoff = now - seconds * 1000;
  return Object.values(jobs)
    .filter((job) => {
      if (!isTerminal(job.status)) return false;
      if (!job.finished_at) return true;
      const finished = Date.parse(job.finished_at);
      return Number.isNaN(finished) ? true : finished >= cutoff;
    })
    .sort((a, b) => byCreatedAt(b, a));
}

/** Everything the console shows: running first, then queued, then just-finished. */
export function selectConsoleJobs(jobs: Record<string, Job>, now?: number): Job[] {
  return [...selectRunning(jobs), ...selectQueued(jobs), ...selectRecentTerminal(jobs, undefined, now)];
}
