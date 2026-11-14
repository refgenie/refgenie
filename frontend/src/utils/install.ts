/** Used when the server does not report its version. */
export const FALLBACK_VERSION = '1.0.0a2';

/**
 * The pip line for the refgenie version this server runs, so it never goes
 * stale. The `>=` floor also lets pip pick a pre-release without `--pre`.
 * `extra` adds an optional set of dependencies, e.g. `dash`. A local build
 * label (`+g1234`) is dropped: pip does not allow one after `>=`.
 */
export function installCommand(version: string | null | undefined, extra?: string): string {
  const known = version && version !== 'unknown' ? version.split('+')[0] : '';
  const floor = known || FALLBACK_VERSION;
  const name = extra ? `refgenie[${extra}]` : 'refgenie';
  return `pip install "${name}>=${floor}"`;
}
