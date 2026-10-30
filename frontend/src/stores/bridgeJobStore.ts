/**
 * Pulls this page handed to a connected local refgenie over the bridge.
 *
 * Deliberately a SEPARATE store from `jobStore`: that one mirrors THIS
 * instance's jobs (fed by its SSE stream, read by its console and its
 * invalidation bus), and a job running on another machine must never enter
 * it. This store only remembers which jobs the page started, on which local
 * instance, and the last record polled for each.
 *
 * Every change is written through to `localStorage`, so a pull keeps
 * reporting after the user navigates away from the button that started it,
 * and after a reload. Polling itself lives in `hooks/useBridgeJobPolling.ts`.
 */

import { create } from 'zustand';
import { isTerminal } from '../services/contracts';
import { recallBridgeJobs, rememberBridgeJobs } from '../services/bridge/persistence';
import type { Job } from '../services/contracts';

/** The oldest entries are dropped past this many, so storage cannot grow forever. */
export const BRIDGE_JOB_CAP = 20;

export interface BridgeJobEntry {
  /** Job id on the local refgenie. */
  id: string;
  /** `/ping` `instance_id` of the local refgenie running it. */
  instanceId: string;
  /** What was pulled, for display before the first poll returns. */
  label: string;
  genomeDigest: string;
  assetGroupName: string;
  assetName: string | null;
  /**
   * The remote server the pull was asked to fetch from, so a "not subscribed"
   * failure can link to the local subscribe form for it. Optional: entries
   * persisted before this field existed lack it.
   */
  serverUrl?: string;
  submittedAt: string;
  /** The last polled record, or null before the first poll. */
  job: Job | null;
  /** The local refgenie no longer knows the job (404: it restarted). */
  lost: boolean;
}

export interface BridgeJobStoreState {
  entries: BridgeJobEntry[];
  track: (entry: Omit<BridgeJobEntry, 'job' | 'lost'> & { job?: Job | null }) => void;
  update: (id: string, job: Job) => void;
  markLost: (id: string) => void;
  dismiss: (id: string) => void;
  clearFinished: (instanceId: string) => void;
  /** Forget everything (an explicit Disconnect). */
  clear: () => void;
}

/** Still worth polling: not lost, and no terminal record yet. */
export function isActiveEntry(entry: BridgeJobEntry): boolean {
  return !entry.lost && (entry.job === null || !isTerminal(entry.job.status));
}

export function entriesFor(
  entries: readonly BridgeJobEntry[],
  instanceId: string | null | undefined,
): BridgeJobEntry[] {
  return instanceId ? entries.filter((entry) => entry.instanceId === instanceId) : [];
}

function isEntry(value: unknown): value is BridgeJobEntry {
  if (typeof value !== 'object' || value === null) return false;
  const v = value as Record<string, unknown>;
  return (
    typeof v.id === 'string' &&
    typeof v.instanceId === 'string' &&
    typeof v.label === 'string' &&
    typeof v.genomeDigest === 'string' &&
    typeof v.assetGroupName === 'string' &&
    typeof v.submittedAt === 'string' &&
    (v.job === null || typeof v.job === 'object')
  );
}

function hydrate(): BridgeJobEntry[] {
  return recallBridgeJobs()
    .filter(isEntry)
    .map((entry) => ({ ...entry, assetName: entry.assetName ?? null, lost: entry.lost === true }));
}

export const useBridgeJobStore = create<BridgeJobStoreState>()((set, get) => {
  const commit = (entries: BridgeJobEntry[]) => {
    set({ entries });
    rememberBridgeJobs(entries);
  };
  const patch = (id: string, change: Partial<BridgeJobEntry>) =>
    commit(get().entries.map((entry) => (entry.id === id ? { ...entry, ...change } : entry)));

  return {
    entries: hydrate(),

    track: (entry) => {
      const rest = get().entries.filter((existing) => existing.id !== entry.id);
      // A duplicate submit returns the existing job id, so replacing by id
      // keeps one card per job rather than one per click.
      commit(
        [{ ...entry, job: entry.job ?? null, lost: false }, ...rest].slice(0, BRIDGE_JOB_CAP),
      );
    },
    update: (id, job) => patch(id, { job }),
    markLost: (id) => patch(id, { lost: true }),
    dismiss: (id) => commit(get().entries.filter((entry) => entry.id !== id)),
    clearFinished: (instanceId) =>
      commit(
        get().entries.filter((entry) => entry.instanceId !== instanceId || isActiveEntry(entry)),
      ),
    clear: () => commit([]),
  };
});
