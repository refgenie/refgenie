/**
 * URLs on this instance's server, derived from the API base so they also work
 * when the SPA and API are on different origins.
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
 * `serverRoot`, always absolute, for giving this server's address to another
 * process (the local refgenie, via a bridge pull's `server_url`).
 */
export function absoluteServerRoot(
  apiBaseUrl: string,
  origin: string = window.location.origin,
): string {
  return new URL(serverRoot(apiBaseUrl) || '/', origin).href.replace(/\/$/, '');
}
