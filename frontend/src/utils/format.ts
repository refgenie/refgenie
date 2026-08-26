/** Formatting helpers shared by every view. */

const UNITS = ['B', 'KB', 'MB', 'GB', 'TB', 'PB'];

/**
 * Human-readable byte size, matching the output of the retired `file_size`
 * Jinja filter. Returns 'NA' for null/undefined so tables stay aligned.
 */
export function formatBytes(bytes: number | null | undefined): string {
  if (bytes === null || bytes === undefined || Number.isNaN(bytes)) return 'NA';
  if (bytes < 1024) return `${bytes} B`;
  let value = bytes;
  let unit = 0;
  while (value >= 1024 && unit < UNITS.length - 1) {
    value /= 1024;
    unit += 1;
  }
  return `${value.toFixed(1)} ${UNITS[unit]}`;
}

/** Truncated digest for display; the full value goes in a `title` attribute. */
export function formatDigest(digest: string | null | undefined, length = 12): string {
  if (!digest) return 'NA';
  return digest.length <= length ? digest : `${digest.slice(0, length)}…`;
}

/** Comma-joined list, or 'NA' when there is nothing to show. */
export function formatList(values: readonly string[] | null | undefined): string {
  if (!values || values.length === 0) return 'NA';
  return values.join(', ');
}
