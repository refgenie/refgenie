/**
 * How each job failure is shown, keyed on `error.code` from
 * `refgenie/server/errors.py`. The older `missing_*` codes are kept as aliases.
 */

export type JobErrorTone = 'warning' | 'info' | 'neutral' | 'error';

/**
 * What the card offers. The ids are resolved to handlers by the component that
 * owns the context (a pull card can re-submit; the jobs history page cannot).
 */
export type JobErrorActionId =
  | 'pull_force'
  | 'pick_server'
  | 'subscribe_server'
  | 'retry'
  | 'view_log'
  | 'fix_inputs'
  | 'build_input'
  | 'init_genome'
  | 'manage_recipes';

export interface JobErrorPresentation {
  tone: JobErrorTone;
  headline: string;
  action?: { id: JobErrorActionId; label: string };
  /** Explanatory copy shown instead of an action. */
  note?: string;
  /** Render `error.detail` as preformatted text (for log tails). */
  preformatted?: boolean;
  /** Show a `<details>` disclosure with `error.detail`. */
  disclosure?: boolean;
}

export const JOB_ERROR_PRESENTATION: Record<string, JobErrorPresentation> = {
  asset_exists: {
    tone: 'warning',
    headline: 'Already downloaded',
    action: { id: 'pull_force', label: 'Pull anyway (overwrite)' },
  },
  already_exists: {
    tone: 'warning',
    headline: 'Already downloaded',
    action: { id: 'pull_force', label: 'Pull anyway (overwrite)' },
  },
  no_archive: {
    tone: 'info',
    headline: 'This server has no archive for that asset',
    action: { id: 'pick_server', label: 'Try another server' },
  },
  server_cannot_serve: {
    tone: 'info',
    headline: 'That server cannot serve this asset',
    action: { id: 'pick_server', label: 'Try another server' },
  },
  no_subscriptions: {
    tone: 'info',
    headline: 'No servers are subscribed',
    action: { id: 'subscribe_server', label: 'Subscribe to a server' },
  },
  server_not_subscribed: {
    tone: 'info',
    headline: 'This refgenie is not subscribed to that server',
    action: { id: 'subscribe_server', label: 'Subscribe to this server' },
  },
  remote_unavailable: {
    tone: 'error',
    headline: 'Could not reach the remote server',
    action: { id: 'retry', label: 'Retry' },
  },
  pull_skipped: {
    tone: 'neutral',
    headline: 'Pull cancelled',
    action: { id: 'retry', label: 'Retry' },
  },
  cancelled: {
    tone: 'neutral',
    headline: 'Cancelled',
    action: { id: 'retry', label: 'Retry' },
  },
  pull_failed: {
    tone: 'error',
    headline: 'Pull failed',
    action: { id: 'retry', label: 'Retry' },
    disclosure: true,
  },
  remote_digest_mismatch: {
    tone: 'error',
    headline: "The server's parent asset differs from your local one",
    note:
      'Pulling would produce an asset whose parents do not match the ones you have. ' +
      'Rebuild or re-pull the parent asset first.',
  },
  // The builder returns None rather than raising, so the detail is the log.
  build_failed: {
    tone: 'error',
    headline: 'Build failed',
    action: { id: 'view_log', label: 'View log' },
    preformatted: true,
  },
  missing_build_input: {
    tone: 'error',
    headline: 'Missing build input',
    action: { id: 'fix_inputs', label: 'Fix inputs' },
  },
  container_engine_missing: {
    tone: 'error',
    headline: 'No container engine',
    note: 'This recipe runs in a container image. Install Docker or Apptainer, then build again.',
  },
  asset_not_found: {
    tone: 'error',
    headline: 'A required input asset is not built yet',
    action: { id: 'build_input', label: 'Build it' },
  },
  asset_group_not_found: {
    tone: 'error',
    headline: 'A required input asset is not built yet',
    action: { id: 'build_input', label: 'Build it' },
  },
  missing_asset: {
    tone: 'error',
    headline: 'A required input asset is not built yet',
    action: { id: 'build_input', label: 'Build it' },
  },
  missing_asset_group: {
    tone: 'error',
    headline: 'A required input asset is not built yet',
    action: { id: 'build_input', label: 'Build it' },
  },
  genome_not_found: {
    tone: 'error',
    headline: 'Genome not found locally',
    action: { id: 'init_genome', label: 'Initialize genome' },
  },
  alias_not_found: {
    tone: 'error',
    headline: 'Genome not found locally',
    action: { id: 'init_genome', label: 'Initialize genome' },
  },
  missing_genome: {
    tone: 'error',
    headline: 'Genome not found locally',
    action: { id: 'init_genome', label: 'Initialize genome' },
  },
  missing_alias: {
    tone: 'error',
    headline: 'Genome not found locally',
    action: { id: 'init_genome', label: 'Initialize genome' },
  },
  recipe_not_found: {
    tone: 'error',
    headline: 'Recipe not found',
    action: { id: 'manage_recipes', label: 'Manage recipes' },
  },
  missing_recipe: {
    tone: 'error',
    headline: 'Recipe not found',
    action: { id: 'manage_recipes', label: 'Manage recipes' },
  },
  asset_class_not_found: {
    tone: 'error',
    headline: 'Asset class not found',
    action: { id: 'manage_recipes', label: 'Manage asset classes' },
  },
};

/** Falls through to a generic retry-plus-disclosure card for unknown codes. */
export function presentJobError(error: {
  code?: string | null;
  message?: string | null;
}): JobErrorPresentation {
  const known = error.code ? JOB_ERROR_PRESENTATION[error.code] : undefined;
  if (known) return known;
  return {
    tone: 'error',
    headline: error.message || 'The job failed',
    action: { id: 'retry', label: 'Retry' },
    disclosure: true,
  };
}
