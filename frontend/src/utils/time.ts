/** Time formatting for job cards and the history table. */

export function formatTimestamp(value: string | null | undefined): string {
  if (!value) return 'NA';
  const parsed = Date.parse(value);
  if (Number.isNaN(parsed)) return value;
  return new Date(parsed).toLocaleString();
}

/** `1m 04s`, or `NA` when the job has not started. */
export function formatElapsed(ms: number): string {
  if (!Number.isFinite(ms) || ms < 0) return 'NA';
  const total = Math.floor(ms / 1000);
  const hours = Math.floor(total / 3600);
  const minutes = Math.floor((total % 3600) / 60);
  const seconds = total % 60;
  const pad = (n: number) => String(n).padStart(2, '0');
  if (hours > 0) return `${hours}h ${pad(minutes)}m`;
  if (minutes > 0) return `${minutes}m ${pad(seconds)}s`;
  return `${seconds}s`;
}

/**
 * Elapsed time between two ISO timestamps. A running job passes `null` as the
 * end and gets the time so far.
 */
export function formatDuration(
  start: string | null | undefined,
  end: string | null | undefined,
  now: number = Date.now(),
): string {
  if (!start) return 'NA';
  const from = Date.parse(start);
  if (Number.isNaN(from)) return 'NA';
  const to = end ? Date.parse(end) : now;
  if (Number.isNaN(to)) return 'NA';
  return formatElapsed(to - from);
}
