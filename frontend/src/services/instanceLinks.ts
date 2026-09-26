/**
 * URLs for endpoints on THIS instance.
 *
 * `client.baseUrl` is `${root_path}/v4` same-origin, or an absolute
 * `https://api.refgenie.org/v4` when the SPA is deployed away from its API
 * (refgenie.org reading api.refgenie.org). Dropping the version suffix gives
 * the server root in both cases, with root_path already folded in by
 * `resolveApiBase`.
 *
 * This is why the pages do not use root-relative literals: `/openapi.json` on
 * the static cross-origin deployment resolves against the SPA origin, where
 * there is no API at all.
 */

/** `https://api.refgenie.org/v4` -> `https://api.refgenie.org`; `/v4` -> ``. */
export function serverRoot(apiBaseUrl: string): string {
  return apiBaseUrl.replace(/\/v\d+$/, '');
}

/** An absolute-or-root-relative URL for a path on this instance's server. */
export function instanceUrl(apiBaseUrl: string, path: string): string {
  const suffix = path.startsWith('/') ? path : `/${path}`;
  return `${serverRoot(apiBaseUrl)}${suffix}`;
}

/**
 * `serverRoot`, but always absolute — for handing this server's address to a
 * DIFFERENT process.
 *
 * The bridge is why this exists. `server_url` in a pull request tells the
 * user's local refgenie which server to fetch the asset from, and that process
 * lives on another origin, where the same-origin deployment's `serverRoot` of
 * `''` names nothing at all. Resolving against the page's own origin is what
 * turns "my API" into an address somebody else can dial.
 */
export function absoluteServerRoot(
  apiBaseUrl: string,
  origin: string = window.location.origin,
): string {
  return new URL(serverRoot(apiBaseUrl) || '/', origin).href.replace(/\/$/, '');
}
