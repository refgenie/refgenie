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
  loc.origin === `http://localhost:${port}` ||
  loc.origin === `http://127.0.0.1:${port}`;

export type PermissionsLike = {
  query: (descriptor: { name: string }) => Promise<{ state: string }>;
};

/**
 * Chromium's names for the permission that gates fetching localhost from a
 * public page. Chrome has renamed it while shipping, so each is tried in turn.
 */
const LOCAL_NETWORK_PERMISSIONS = [
  'loopback-network',
  'local-network-access',
  'local-network',
];

/**
 * Whether probing localhost now would show the browser's permission prompt.
 * `unsupported` means this browser does not gate the request by permission,
 * so a probe is silent.
 */
export const localNetworkPermission = async (
  permissions: PermissionsLike | undefined = typeof navigator !== 'undefined'
    ? (navigator.permissions as unknown as PermissionsLike | undefined)
    : undefined,
): Promise<'granted' | 'prompt' | 'denied' | 'unsupported'> => {
  if (!permissions?.query) return 'unsupported';
  for (const name of LOCAL_NETWORK_PERMISSIONS) {
    try {
      const { state } = await permissions.query({ name });
      if (state === 'granted' || state === 'prompt' || state === 'denied') return state;
    } catch {
      // Unknown permission name in this browser; try the next one.
    }
  }
  return 'unsupported';
};
