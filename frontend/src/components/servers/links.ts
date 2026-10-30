/**
 * Where a user goes to subscribe to a server.
 *
 * Subscribing is a same-origin action on the local dash: a page reached over
 * the bridge may read the subscription list but never change it. So every
 * "subscribe to this server" affordance, local or bridged, is a link to the
 * local dash's Servers page with the form already filled in, never a request.
 */

/**
 * The server a fresh `refgenie init` subscribes to. Mirrors
 * `DEFAULT_SERVER_URL` in `refgenie/const.py`; used to prefill the form when a
 * pull failed without naming a server.
 */
export const DEFAULT_SERVER_URL = 'https://api.refgenie.org';

/** The query parameter that opens the Servers page's form prefilled. */
export const SUBSCRIBE_PARAM = 'subscribe';

/**
 * `/servers?subscribe=<url>`, relative to the app's own route base, or the bare
 * `/servers` without a URL.
 */
export function serversPath(subscribeUrl?: string): string {
  return subscribeUrl
    ? `/servers?${SUBSCRIBE_PARAM}=${encodeURIComponent(subscribeUrl)}`
    : '/servers';
}

/** The same page on a local dash reached over the bridge (`http://localhost:8080`). */
export function localServersUrl(localBaseUrl: string, subscribeUrl?: string): string {
  return `${localBaseUrl.replace(/\/$/, '')}${serversPath(subscribeUrl)}`;
}
