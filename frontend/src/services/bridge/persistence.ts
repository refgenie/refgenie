/**
 * localStorage for the localhost bridge: the remembered connection (which allows
 * a silent re-probe on later visits; the first probe is always user-initiated),
 * the dismissed Safari card, and in-flight pulls. Do not rename the keys or the
 * stored shape, or returning visitors lose them. Every accessor swallows storage
 * errors (private modes, disabled storage).
 */

const CONNECTION_KEY = 'refgenie.localBridge';
/** Pulls started on the local refgenie, so progress survives a reload. */
const JOBS_KEY = 'refgenie.localBridge.jobs';

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

/** Read the remembered bridge jobs, or `[]`. `stores/bridgeJobStore.ts` validates each entry. */
export const recallBridgeJobs = (
  storage: StorageLike | null = defaultStorage(),
): unknown[] => {
  try {
    const raw = storage?.getItem(JOBS_KEY);
    if (!raw) return [];
    const parsed: unknown = JSON.parse(raw);
    return Array.isArray(parsed) ? parsed : [];
  } catch {
    return [];
  }
};

export const rememberBridgeJobs = (
  entries: readonly unknown[],
  storage: StorageLike | null = defaultStorage(),
): void => {
  try {
    if (entries.length === 0) storage?.removeItem(JOBS_KEY);
    else storage?.setItem(JOBS_KEY, JSON.stringify(entries));
  } catch {
    // Storage unavailable: progress just will not survive a reload.
  }
};
