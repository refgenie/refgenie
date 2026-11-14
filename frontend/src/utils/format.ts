/** Formatting helpers shared by every view. */

const UNITS = ['B', 'KB', 'MB', 'GB', 'TB', 'PB'];

/** Human-readable byte size, or 'NA' for a missing value. */
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

/** Base pairs, short: 16,569 bp; 46.7 Mbp; 3.1 Gbp. */
export function formatBp(bases: number): string {
  if (bases >= 1e9) return `${(bases / 1e9).toFixed(1)} Gbp`;
  if (bases >= 1e6) return `${(bases / 1e6).toFixed(1)} Mbp`;
  if (bases >= 1e4) return `${(bases / 1e3).toFixed(1)} kbp`;
  return `${bases.toLocaleString()} bp`;
}
