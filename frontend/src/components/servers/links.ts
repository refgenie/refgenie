/**
 * Links to the local dash's Servers page with the subscribe form prefilled.
 * Subscribing is same-origin only, so bridged pages link rather than request.
 */

/** The server a fresh `refgenie init` subscribes to (`DEFAULT_SERVER_URL` in `refgenie/const.py`). */
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
