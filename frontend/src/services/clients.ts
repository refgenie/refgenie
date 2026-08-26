/**
 * Client construction and the React context that carries the resolved one.
 *
 * There is no module-level singleton and no module-level API base constant:
 * that is the specific defect this layer exists to avoid. A second backend is
 * `createClient(otherUrl)`.
 */

import { createContext } from 'react';
import { ApiClient } from './http';
import type { UiConfig } from '../types/ui';

export function createClient(baseUrl: string): ApiClient {
  return new ApiClient({ baseUrl });
}

/**
 * Resolve the effective API base from the bootstrap config: a relative
 * `api_base` is prefixed with `root_path` so sub-path deployments work.
 */
export function resolveApiBase(config: Pick<UiConfig, 'api_base' | 'root_path'>): string {
  const base = config.api_base;
  if (/^https?:\/\//i.test(base)) return base;
  const root = (config.root_path ?? '').replace(/\/+$/, '');
  const path = base.startsWith('/') ? base : `/${base}`;
  return `${root}${path}`;
}

/**
 * Base for the local-only command/remote surface (`/v1`). It is always
 * same-origin: unlike `api_base` it is never cross-origin, because a local
 * dash is the only thing that serves it.
 */
export function resolveLocalApiBase(config: Pick<UiConfig, 'root_path'>): string {
  return `${(config.root_path ?? '').replace(/\/+$/, '')}/v1`;
}

export interface AppContextValue {
  config: UiConfig;
  /** Shared read API (`/v4`). */
  client: ApiClient;
  /** Local-only surface (`/v1`): remote browse, and later actions and jobs. */
  localClient: ApiClient;
}

/** Undefined only outside the provider, which is a programming error. */
export const AppContext = createContext<AppContextValue | undefined>(undefined);
