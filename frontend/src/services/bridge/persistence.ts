/**
 * Remembered-connection persistence for the localhost bridge.
 *
 * `localStorage` key `refgenie.localBridge` stores `{port, instanceId,
 * connectedAt}` after a successful connect; its presence licenses a silent
 * re-probe on later visits (D3 — the first probe is always user-initiated).
 *
 * Both key names and the stored JSON shape are frozen: refgenie.org is served
 * today by the refgenie-ui SPA, and this app will replace it on the same
 * origin. A returning visitor's remembered connection and dismissed Safari card
 * carry across that swap only as long as nothing here is renamed.
 *
 * Storage is injectable for tests; every accessor is exception-safe because
 * localStorage can throw (privacy modes, disabled storage).
 */

const CONNECTION_KEY = 'refgenie.localBridge';
const SAFARI_DISMISS_KEY = 'refgenie.localBridge.safariDismissed';

/** The `refgenie dash` default port, and the vite dev proxy target. */
export const DEFAULT_LOCAL_PORT = 8080;

export type RememberedConnection = {
  port: number;
  instanceId: string;
  connectedAt: string;
};

type StorageLike = Pick<Storage, 'getItem' | 'setItem' | 'removeItem'>;

const defaultStorage = (): StorageLike | null => {
  try {
    return typeof window !== 'undefined' ? window.localStorage : null;
  } catch {
    return null;
  }
};

export const remember = (
  connection: RememberedConnection,
  storage: StorageLike | null = defaultStorage(),
): void => {
  try {
    storage?.setItem(CONNECTION_KEY, JSON.stringify(connection));
  } catch {
    // Storage unavailable: the connection just will not be remembered.
  }
};

export const recall = (
  storage: StorageLike | null = defaultStorage(),
): RememberedConnection | null => {
  try {
    const raw = storage?.getItem(CONNECTION_KEY);
    if (!raw) return null;
    const parsed: unknown = JSON.parse(raw);
    if (
      typeof parsed === 'object' &&
      parsed !== null &&
      typeof (parsed as RememberedConnection).port === 'number' &&
      typeof (parsed as RememberedConnection).instanceId === 'string'
    ) {
      return parsed as RememberedConnection;
    }
    return null;
  } catch {
    return null;
  }
};

export const forget = (storage: StorageLike | null = defaultStorage()): void => {
  try {
    storage?.removeItem(CONNECTION_KEY);
  } catch {
    // Nothing to clear.
  }
};

export const isSafariCardDismissed = (
  storage: StorageLike | null = defaultStorage(),
): boolean => {
  try {
    return storage?.getItem(SAFARI_DISMISS_KEY) === 'true';
  } catch {
    return false;
  }
};

export const dismissSafariCard = (
  storage: StorageLike | null = defaultStorage(),
): void => {
  try {
    storage?.setItem(SAFARI_DISMISS_KEY, 'true');
  } catch {
    // Best effort.
  }
};
