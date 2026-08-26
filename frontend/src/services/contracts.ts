/**
 * The single file to edit when a sibling backend plan's shape changes.
 *
 * Everything the management surface needs from the local-only `/v1` API lives
 * here: the endpoint paths, the job record, the SSE frame shapes, and the
 * action request/response bodies. Components and hooks import from here and
 * never hard-code a path.
 *
 * Paths are RELATIVE to the local client's base (`/v1`), because that base is
 * resolved from `/service-info` at startup — see `services/clients.ts`. There
 * is deliberately no `/v1` literal in this file's path constants.
 */

export { ACTION_HEADER, ACTION_HEADER_VALUE } from './http';

// === Endpoint paths (relative to the `/v1` local base) ===

export const ACTION_PATHS = {
  pull: '/actions/pull',
  build: '/actions/build',
  buildPreflight: '/actions/build/preflight',
  genomes: '/actions/genomes',
  aliases: '/actions/aliases',
  subscriptions: '/actions/subscriptions',
  assetDefault: '/actions/assets/default',
  asset: (digest: string) => `/actions/assets/${encodeURIComponent(digest)}`,
  genome: (ref: string) => `/actions/genomes/${encodeURIComponent(ref)}`,
  alias: (name: string) => `/actions/aliases/${encodeURIComponent(name)}`,
} as const;

export const JOB_PATHS = {
  list: '/jobs',
  events: '/jobs/events',
  job: (id: string) => `/jobs/${encodeURIComponent(id)}`,
  log: (id: string) => `/jobs/${encodeURIComponent(id)}/log`,
  cancel: (id: string) => `/jobs/${encodeURIComponent(id)}/cancel`,
} as const;

/**
 * The `action` label passed to `ApiClient.mutate`. It is not sent as the header
 * value (the server ignores that), but it names the operation at every call
 * site and keeps the set of mutations enumerable.
 */
export type ActionVerb =
  | 'pull'
  | 'build'
  | 'build.preflight'
  | 'genome.init'
  | 'asset.delete'
  | 'genome.delete'
  | 'alias.set'
  | 'alias.remove'
  | 'subscribe'
  | 'unsubscribe'
  | 'asset.set_default'
  | 'job.cancel';

// === Jobs ===

export type JobKind = 'pull' | 'build' | 'genome_init';

/**
 * `pending` and `completed` are not statuses. `queued` is load-bearing: builds
 * run in a 1-slot executor and pulls in a 2-slot one, so a submitted job
 * routinely waits and the UI must say so with a position rather than a spinner.
 */
export type JobStatus = 'queued' | 'running' | 'succeeded' | 'failed' | 'cancelled';

export const TERMINAL_STATUSES: readonly JobStatus[] = ['succeeded', 'failed', 'cancelled'];

export function isTerminal(status: JobStatus): boolean {
  return TERMINAL_STATUSES.includes(status);
}

/** Structured job identity; the basis of the duplicate/active-job lookup. */
export interface JobTarget {
  genome_digest: string | null;
  genome_name: string | null;
  asset_group_name: string;
  asset_name: string | null;
}

export interface JobProgressInfo {
  phase: string;
  /** `null` for every build and for any pull the manager did not instrument. */
  percent: number | null;
  message: string | null;
  bytes_done: number | null;
  /** Legitimately `null`: legacy staged rows carry no `tarball_size`. */
  bytes_total: number | null;
}

export interface JobResult {
  asset_digest: string;
  registry_path: string;
}

export interface JobError {
  code: string;
  message: string;
  detail: string | null;
  field: string | null;
}

export interface Job {
  id: string;
  kind: JobKind;
  status: JobStatus;
  /** Position within this job's OWN executor queue; `null` once it starts. */
  queue_position: number | null;
  label: string;
  target: JobTarget;
  progress: JobProgressInfo | null;
  created_at: string;
  started_at: string | null;
  finished_at: string | null;
  result: JobResult | null;
  error: JobError | null;
  log_lines: number;
  /** False when a running build is past the point of cooperative cancellation. */
  cancellable?: boolean;
  /**
   * A build whose asset already existed: it succeeded without doing any work
   * (`AssetBuilder.build` restores the flag and returns the existing asset).
   */
  skipped?: boolean;
  /** Resolved build commands, when the record carries them (builds only). */
  build_commands?: string[] | null;
}

/** The 202 body of every job-producing action. */
export interface JobRef {
  job_id: string;
  kind: JobKind;
  status: JobStatus;
  created_at: string;
  /**
   * True when an identical submission was already in flight. The manager
   * coalesces rather than rejecting, so there is no 409 branch anywhere.
   */
  duplicate: boolean;
  links?: { self?: string; events?: string; cancel?: string };
}

export interface JobLogPage {
  lines: string[];
  next_offset: number;
  truncated: boolean;
}

// === SSE frames ===

export type JobEventType = 'status' | 'progress' | 'log' | 'done' | 'truncated' | 'heartbeat';

interface JobEventBase {
  /** Manager-global and monotonic across ALL jobs, not per-job. */
  seq: number;
  job_id?: string;
}

export interface JobStatusEvent extends JobEventBase {
  type: 'status';
  job_id: string;
  status: JobStatus;
  queue_position?: number | null;
  job?: Job;
}

export interface JobProgressEvent extends JobEventBase, Partial<JobProgressInfo> {
  type: 'progress';
  job_id: string;
}

export interface JobLogEvent extends JobEventBase {
  type: 'log';
  job_id: string;
  line?: string;
  lines?: string[];
  /** `pipeline` is tee'd subprocess output; anything else is refgenie's logger. */
  source?: 'refgenie' | 'pipeline' | string;
}

export interface JobDoneEvent extends JobEventBase {
  type: 'done';
  job_id: string;
  /** The full record. Some servers spread it at the top level instead. */
  job?: Job;
}

export interface JobTruncatedEvent extends JobEventBase {
  type: 'truncated';
  job_id: string;
  dropped?: number;
}

export interface JobHeartbeatEvent extends JobEventBase {
  type: 'heartbeat';
  ts?: string;
}

export type JobEvent =
  | JobStatusEvent
  | JobProgressEvent
  | JobLogEvent
  | JobDoneEvent
  | JobTruncatedEvent
  | JobHeartbeatEvent;

export const JOB_EVENT_TYPES: readonly JobEventType[] = [
  'status',
  'progress',
  'log',
  'done',
  'truncated',
  'heartbeat',
];

/** Named heartbeat interval the server promises; the watchdog is 3x this. */
export const HEARTBEAT_INTERVAL_MS = 15_000;

// === Action request / response bodies ===

/*
 * Every request model below is `extra="forbid"` server-side
 * (`refgenie/server/local/schemas.py`), so a stale or invented field name is a
 * loud 422 rather than a silent no-op. Field names here are transcribed from
 * that module and must match it exactly.
 */

export interface PullRequest {
  asset_group: string;
  /** EXACTLY ONE of `genome` (an alias) or `genome_digest` may be set. */
  genome?: string | null;
  genome_digest?: string | null;
  asset?: string | null;
  server_url?: string | null;
  /**
   * ALWAYS sent explicitly. A tri-state `force` means "ask the user", and
   * `force: null` server-side reaches a `Confirm.ask` on stdin that hangs the
   * worker forever, so the client never omits it.
   */
  force: boolean;
}

/** The recipe-driven inputs, NESTED under `params` on a build request. */
export interface BuildParamsRequest {
  assets?: Record<string, string> | null;
  params?: Record<string, string | number | boolean> | null;
  files?: Record<string, string> | null;
}

export interface BuildRequest {
  recipe: string;
  /** An ALIAS. `preflight_build` resolves it via `alias.resolve`, which does
   *  not accept a digest. */
  genome: string;
  asset_group: string;
  asset?: string | null;
  recipe_version?: string | null;
  description?: string | null;
  stage?: boolean;
  pull_parents?: boolean;
  params?: BuildParamsRequest | null;
}

export interface PreflightFieldError {
  /** Always present. Dotted, matching the request body's own field paths. */
  field: string;
  code: string;
  message: string;
}

/**
 * `resolved` reports what the build WOULD use. There is no `registry_path`:
 * the caller composes it from the genome, the asset group and `asset_name`.
 */
export interface PreflightResolved {
  genome_digest?: string;
  recipe?: string;
  recipe_version?: string;
  input_assets?: Record<string, string | null>;
  asset_name?: string;
  [key: string]: unknown;
}

export interface PreflightResponse {
  /** Always a 200: a preflight that found problems has succeeded at its job. */
  ok: boolean;
  errors: PreflightFieldError[];
  resolved: PreflightResolved;
}

export interface GenomeInitRequest {
  /** A server-local path or a URL. */
  fasta: string;
  aliases: string[];
  description?: string | null;
  species?: string | null;
  build_fasta_asset: boolean;
}

export interface AliasRequest {
  alias: string;
  /** Required: the digest-less facade variant does synchronous network
   *  lookups against every subscription, which a request/response endpoint
   *  must not do. */
  genome_digest: string;
}

export interface SubscribeRequest {
  /** Always the list form, and never empty (`min_length=1`). */
  server_urls: string[];
  reset?: boolean;
}

export interface UnsubscribeRequest {
  server_urls: string[];
}

export interface SetDefaultAssetRequest {
  genome_digest: string;
  asset_group: string;
  asset: string;
}

/** The 200 body of every synchronous action. */
export interface ActionResult {
  ok: true;
  message: string;
  data: Record<string, unknown> | null;
}
