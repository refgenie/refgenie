/** Browser checks for the localhost bridge; inputs are injectable for tests. */

export type NavigatorLike = { vendor?: string; userAgent?: string };
export type LocationLike = { origin: string };

/**
 * Safari or any iOS browser. WebKit blocks HTTPS pages from fetching
 * http://localhost (WebKit bug 171934), so the probe is skipped there.
 */
export const isWebKitOnly = (
  nav: NavigatorLike = typeof navigator !== 'undefined' ? navigator : {},
): boolean => {
  const userAgent = nav.userAgent ?? '';
  const chromiumMarkers = /Chrome|Chromium|CriOS|Edg\/|OPR\//;
  return nav.vendor === 'Apple Computer, Inc.' && !chromiumMarkers.test(userAgent);
};

/**
 * Whether this page is the local dash itself. Backs up the `mode` check for a
 * dash whose `/service-info` failed, so it does not offer to connect to itself.
 */
export const isSelfLocal = (
  port: number,
  loc: LocationLike = typeof window !== 'undefined' ? window.location : { origin: '' },
): boolean =>
  loc.origin === `http://localhost:${port}` || loc.origin === `http://127.0.0.1:${port}`;
