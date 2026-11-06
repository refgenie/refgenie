/**
 * Wire contract for the local-only `/v1` API: endpoint paths, job records, SSE
 * frames, and action bodies. Paths are relative to the `/v1` base, which is
 * resolved from `/service-info` at startup (see `services/clients.ts`).
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
  dataChannels: '/actions/data_channels',
  dataChannel: (name: string) => `/actions/data_channels/${encodeURIComponent(name)}`,
  dataChannelSync: (name: string) => `/actions/data_channels/${encodeURIComponent(name)}/sync`,
  assetDefault: '/actions/assets/default',
  asset: (digest: string) => `/actions/assets/${encodeURIComponent(digest)}`,
  genome: (digest: string) => `/actions/genomes/${encodeURIComponent(digest)}`,
  alias: (name: string) => `/actions/aliases/${encodeURIComponent(name)}`,
} as const;

export const JOB_PATHS = {
  list: '/jobs',
  events: '/jobs/events',
  job: (id: string) => `/jobs/${encodeURIComponent(id)}`,
  log: (id: string) => `/jobs/${encodeURIComponent(id)}/log`,
  cancel: (id: string) => `/jobs/${encodeURIComponent(id)}/cancel`,
} as const;

/** Names of the mutations passed as `action` to `ApiClient.mutate`. */
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
  | 'data_channel.add'
  | 'data_channel.remove'
  | 'data_channel.sync'
  | 'asset.set_default'
  | 'job.cancel';

// === Jobs ===

export type JobKind = 'pull' | 'build' | 'genome_init';

/**
 * Builds run one at a time and pulls two at a time, so a new job often waits;
 * the UI shows its `queued` position rather than a spinner.
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
  /** `null` for file-mode (non-archive) staged rows, which have no size. */
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
  /** A build whose asset already existed, so it succeeded without doing work. */
  skipped?: boolean;
}

/** The 202 body of every job-producing action. */
export interface JobRef {
  job_id: string;
  kind: JobKind;
  status: JobStatus;
  created_at: string;
  /** True when an identical job was already in flight; the server reuses it instead of returning 409. */
  duplicate: boolean;
  links?: { self?: string; events?: string; cancel?: string };
}

export interface JobLogPage {
  lines: string[];
  next_offset: number;
  /** More lines exist after this page; fetch again from `next_offset`. */
  has_more: boolean;
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

// === Action request / response bodies ===

/*
 * Field names must match `refgenie/server/jobs/schemas.py` (pull, build, genome
 * init) and `refgenie/server/local/schemas.py` (the rest) exactly; the server
 * rejects unknown fields with a 422.
 */

export interface PullRequest {
  asset_group: string;
  /** Set exactly one of `genome` (an alias) or `genome_digest`. */
  genome?: string | null;
  genome_digest?: string | null;
  asset?: string | null;
  server_url?: string | null;
  /** Always sent: `force: null` makes the server prompt on stdin, which hangs the worker. */
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
  /** An alias, not a digest; preflight resolves it by alias only. */
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

/** What the build would use. The caller builds the registry path from these fields. */
export interface PreflightResolved {
  genome_digest?: string;
  recipe?: string;
  recipe_version?: string;
  input_assets?: Record<string, string | null>;
  asset_name?: string;
  [key: string]: unknown;
}

export interface PreflightResponse {
  /** Preflight returns 200 even when it finds problems; check `ok`. */
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
  /** Required so the server does not look the alias up on every subscribed server. */
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

export interface DataChannelAddRequest {
  name: string;
  /** http(s) only; the server rejects local paths. */
  index_address: string;
  description?: string | null;
  /** Register the channel's recipes and asset classes right away (server default: true). */
  sync?: boolean;
}

/** `data.sync` of an add or sync action: what the sync did, item by item. */
export interface DataChannelSyncReport {
  channel: string;
  asset_classes_added: number;
  asset_classes_skipped: number;
  asset_classes_failed: number;
  recipes_added: number;
  recipes_skipped: number;
  recipes_failed: number;
  errors: string[];
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
