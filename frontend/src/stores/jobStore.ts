/**
 * Live job state. A zustand store, not context, because log lines arrive many
 * times a second and selectors let only the log tail re-render.
 */

import { create } from 'zustand';
import { appCache } from '../app/resourceCache';
import { staleAfter } from '../services/queryKeys';
import { isTerminal } from '../services/contracts';
import type {
  Job,
  JobEvent,
  JobKind,
  JobProgressInfo,
  JobRef,
  JobTarget,
} from '../services/contracts';
import type { ToastApi } from '../components/common/toastContext';

/** Per-job log ring buffer cap. The full log comes from `/v1/jobs/{id}/log`. */
export const LOG_RING_CAP = 500;

/** How long a finished job stays on the console so a fast pull can be read. */
const TERMINAL_LINGER_SECONDS = 60;

export type ConnectionState = 'connecting' | 'open' | 'polling';

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
  focusJob: (jobId: string | null) => void;
  activeJobFor: (targetKey: string) => Job | null;
  clearTerminal: () => void;
  reset: () => void;
}

/** Identity for duplicate/active jobs. Prefers the genome digest so alias and digest submits match. */
export function targetKey(kind: JobKind, target: JobTarget | null | undefined): string {
  const genome = target?.genome_digest ?? target?.genome_name ?? '?';
  const group = target?.asset_group_name ?? '?';
  const name = target?.asset_name ?? '*';
  return `${kind}:${genome}/${group}:${name}`;
}

function jobTargetKey(job: Job): string {
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

  registerQueued: (job) =>
    set((state) => ({
      jobs: { ...state.jobs, [job.id]: { ...state.jobs[job.id], ...job } },
      focusedJobId: job.id,
    })),

  applyEvent: (event) => {
    const seq = typeof event.seq === 'number' ? event.seq : 0;
    // Ignore replayed events at or below the cursor (reconnects may overlap).
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
        for (const family of staleAfter[finished.kind] ?? []) appCache.invalidate([family]);
      }
    }
  },

  appendLog: (jobId, lines, source = 'refgenie') =>
    set((state) => ({ logs: { ...state.logs, [jobId]: pushLines(state.logs[jobId], lines, source) } })),

  setConnection: (connection) => set({ connection }),

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

/** Queued jobs, grouped by kind, then by queue position (each kind has its own queue). */
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

/**
 * Puts a just-submitted job on the console before its first SSE frame, or
 * focuses the card already running when the server reports a duplicate.
 */
export function trackJob(
  ref: JobRef,
  { noun, label, target }: { noun: string; label: string; target: JobTarget },
  toast: ToastApi,
): void {
  const store = useJobStore.getState();
  if (ref.duplicate) {
    store.focusJob(ref.job_id);
    toast.info(`That ${noun} is already running.`);
    return;
  }
  store.registerQueued({
    id: ref.job_id,
    kind: ref.kind,
    status: ref.status,
    queue_position: null,
    label,
    target,
    progress: null,
    created_at: ref.created_at,
    started_at: null,
    finished_at: null,
    result: null,
    error: null,
    log_lines: 0,
    cancellable: true,
  });
  toast.success(
    `${noun[0].toUpperCase()}${noun.slice(1)} queued. Progress is in the job console.`,
  );
}
