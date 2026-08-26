/**
 * The single job-event transport. Mounted once, by `AppLayout`, gated on
 * `capabilities.jobs`. Nothing else in the tree opens a connection.
 *
 * Two transport facts drive the whole design:
 *  - `EventSource` cannot send request headers, so the stream is a plain GET
 *    and `X-Refgenie-Action` lives only on mutations.
 *  - Browsers cap ~6 HTTP/1.1 connections per origin, so there is exactly ONE
 *    multiplexed stream for every job, never one per job.
 *
 * And one failure fact: a silently stalled proxy never fires `onerror`, so a
 * named `heartbeat` event (not an SSE comment, which `EventSource` never
 * dispatches to JS) feeds a watchdog that reconnects on silence.
 */

import { useCallback, useEffect, useRef } from 'react';
import { useLocalApiClient } from './useApiClient';
import { useJobStore } from '../stores/jobStore';
import { jobEventsUrl, listJobs } from '../services/jobs';
import { HEARTBEAT_INTERVAL_MS, JOB_EVENT_TYPES } from '../services/contracts';
import type { JobEvent, JobEventType } from '../services/contracts';
import type { JobsFetcher } from '../services/jobs';
import type { ConnectionState } from '../stores/jobStore';

/**
 * The two methods plus two slots a fake needs to implement exactly. Narrower
 * than the DOM `EventSource` on purpose: anything wider invites a test double
 * that diverges from what the hook actually uses.
 */
export interface EventSourceLike {
  addEventListener(type: string, listener: (event: MessageEvent) => void): void;
  close(): void;
  readyState?: number;
  onopen: ((event: Event) => void) | null;
  onerror: ((event: Event) => void) | null;
}

export interface UseJobEventsOptions {
  enabled?: boolean;
  /** Test seam. Defaults to the browser `EventSource`. */
  eventSourceFactory?: (url: string) => EventSourceLike;
  /** Test seam. Defaults to `listJobs` bound to the local client. */
  fetchJobs?: JobsFetcher;
  backoffMs?: number[];
  /** Consecutive failures before giving up on SSE and polling instead. */
  failoverAfter?: number;
  pollIntervalMs?: { active: number; idle: number };
  heartbeatMs?: number;
  /** How often to try SSE again while polling. */
  sseRetryMs?: number;
}

export interface UseJobEventsResult {
  connection: ConnectionState;
  reconnectNow: () => void;
}

const DEFAULT_BACKOFF = [1000, 2000, 4000, 8000, 15000];
const DEFAULT_FAILOVER_AFTER = 3;
const DEFAULT_POLL = { active: 2000, idle: 15000 };
const DEFAULT_SSE_RETRY_MS = 60_000;
/** Terminal sweep size per poll cycle, so a job that finished between polls lands. */
const TERMINAL_SWEEP_LIMIT = 20;
const ACTIVE_PAGE_LIMIT = 100;
const HYDRATE_TERMINAL_LIMIT = 50;

type Timer = ReturnType<typeof setTimeout>;

export function useJobEvents(options: UseJobEventsOptions = {}): UseJobEventsResult {
  const client = useLocalApiClient();
  const connection = useJobStore((state) => state.connection);
  const enabled = options.enabled ?? true;

  const optionsRef = useRef(options);
  optionsRef.current = options;
  const reconnectRef = useRef<() => void>(() => {});

  useEffect(() => {
    if (!enabled) return;

    const store = () => useJobStore.getState();
    const opts = () => optionsRef.current;
    const fetchJobs: JobsFetcher = (params, init) =>
      (opts().fetchJobs ?? ((p, i) => listJobs(client, p, i)))(params, init);

    let disposed = false;
    let source: EventSourceLike | null = null;
    let failures = 0;
    let polling = false;
    let reconnectTimer: Timer | undefined;
    let watchdogTimer: Timer | undefined;
    let pollTimer: Timer | undefined;
    let sseRetryTimer: Timer | undefined;

    const clear = (timer: Timer | undefined) => {
      if (timer !== undefined) clearTimeout(timer);
    };

    const closeSource = () => {
      clear(watchdogTimer);
      watchdogTimer = undefined;
      if (source) {
        source.onopen = null;
        source.onerror = null;
        source.close();
        source = null;
      }
    };

    const stopPolling = () => {
      polling = false;
      clear(pollTimer);
      pollTimer = undefined;
      clear(sseRetryTimer);
      sseRetryTimer = undefined;
    };

    const stopEverything = () => {
      closeSource();
      stopPolling();
      clear(reconnectTimer);
      reconnectTimer = undefined;
    };

    const hasActiveJob = () =>
      Object.values(store().jobs).some(
        (job) => job.status === 'queued' || job.status === 'running',
      );

    const hydrate = async () => {
      try {
        const active = await fetchJobs({ status: 'active', limit: ACTIVE_PAGE_LIMIT });
        if (disposed) return;
        store().upsertMany(active.items ?? []);
        const terminal = await fetchJobs({
          status: 'terminal',
          limit: HYDRATE_TERMINAL_LIMIT,
        });
        if (disposed) return;
        store().upsertMany(terminal.items ?? []);
      } catch {
        // A failed hydrate is not fatal: the stream still carries everything
        // from the cursor forward.
      }
    };

    const pollOnce = async () => {
      try {
        const active = await fetchJobs({ status: 'active', limit: ACTIVE_PAGE_LIMIT });
        if (disposed) return;
        store().upsertMany(active.items ?? []);
        const terminal = await fetchJobs({ status: 'terminal', limit: TERMINAL_SWEEP_LIMIT });
        if (disposed) return;
        store().upsertMany(terminal.items ?? []);
      } catch {
        // Keep polling: the next cycle may succeed.
      }
    };

    const schedulePoll = () => {
      clear(pollTimer);
      const intervals = opts().pollIntervalMs ?? DEFAULT_POLL;
      const wait = hasActiveJob() ? intervals.active : intervals.idle;
      pollTimer = setTimeout(() => {
        void pollOnce().then(() => {
          if (!disposed && polling) schedulePoll();
        });
      }, wait);
    };

    const scheduleSseRetry = () => {
      clear(sseRetryTimer);
      sseRetryTimer = setTimeout(() => {
        if (disposed || !polling) return;
        openSource({ quiet: true });
      }, opts().sseRetryMs ?? DEFAULT_SSE_RETRY_MS);
    };

    const startPolling = () => {
      if (polling) {
        scheduleSseRetry();
        return;
      }
      polling = true;
      store().setConnection('polling');
      void pollOnce().then(() => {
        if (!disposed && polling) schedulePoll();
      });
      scheduleSseRetry();
    };

    const resetWatchdog = () => {
      clear(watchdogTimer);
      const heartbeat = opts().heartbeatMs ?? HEARTBEAT_INTERVAL_MS;
      watchdogTimer = setTimeout(() => {
        if (disposed) return;
        // Silence, not an error: browsers do not fire `onerror` for a proxy
        // that keeps the socket open and sends nothing.
        handleFailure();
      }, heartbeat * 3);
    };

    const handleFrame = (type: JobEventType, event: MessageEvent) => {
      resetWatchdog();
      let payload: unknown;
      try {
        payload = JSON.parse(String(event.data));
      } catch {
        // One bad line must never kill the stream.
        console.warn('jobs: skipped a malformed SSE frame', event.data);
        return;
      }
      if (!payload || typeof payload !== 'object') return;
      store().applyEvent({ ...(payload as object), type } as JobEvent);
    };

    const handleFailure = () => {
      if (disposed) return;
      closeSource();
      failures += 1;
      const failoverAfter = opts().failoverAfter ?? DEFAULT_FAILOVER_AFTER;
      if (failures >= failoverAfter) {
        startPolling();
        return;
      }
      store().setConnection('offline');
      const backoff = opts().backoffMs ?? DEFAULT_BACKOFF;
      const wait = backoff[Math.min(failures - 1, backoff.length - 1)];
      clear(reconnectTimer);
      reconnectTimer = setTimeout(() => {
        if (!disposed) openSource();
      }, wait);
    };

    function openSource(args: { quiet?: boolean } = {}) {
      if (disposed) return;
      closeSource();
      const factory =
        opts().eventSourceFactory ??
        (typeof EventSource === 'undefined'
          ? undefined
          : (url: string) => new EventSource(url) as unknown as EventSourceLike);

      if (!factory) {
        // No SSE in this environment at all: poll rather than pretend.
        startPolling();
        return;
      }

      if (!args.quiet) store().setConnection('connecting');
      const url = jobEventsUrl(client, store().lastSeq);
      const next = factory(url);
      source = next;

      next.onopen = () => {
        if (disposed || source !== next) return;
        failures = 0;
        stopPolling();
        store().setConnection('open');
        resetWatchdog();
      };
      next.onerror = () => {
        if (disposed || source !== next) return;
        handleFailure();
      };
      for (const type of JOB_EVENT_TYPES) {
        next.addEventListener(type, (event) => {
          if (disposed || source !== next) return;
          handleFrame(type, event);
        });
      }
      resetWatchdog();
    }

    const onVisibilityChange = () => {
      if (typeof document === 'undefined') return;
      if (document.visibilityState === 'hidden') {
        // A background tab must not hold a connection open indefinitely.
        stopEverything();
        store().setConnection('offline');
      } else {
        failures = 0;
        void hydrate().then(() => {
          if (!disposed) openSource();
        });
      }
    };

    reconnectRef.current = () => {
      failures = 0;
      stopPolling();
      clear(reconnectTimer);
      reconnectTimer = undefined;
      openSource();
    };

    void hydrate().then(() => {
      if (!disposed) openSource();
    });

    if (typeof document !== 'undefined') {
      document.addEventListener('visibilitychange', onVisibilityChange);
    }

    return () => {
      disposed = true;
      stopEverything();
      if (typeof document !== 'undefined') {
        document.removeEventListener('visibilitychange', onVisibilityChange);
      }
      reconnectRef.current = () => {};
    };
  }, [enabled, client]);

  const reconnectNow = useCallback(() => reconnectRef.current(), []);

  return { connection, reconnectNow };
}
