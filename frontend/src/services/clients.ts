/** API client construction and the React context that carries the clients. */

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

/** Base for the local-only `/v1` API, always same-origin. */
export function resolveLocalApiBase(config: Pick<UiConfig, 'root_path'>): string {
  return `${(config.root_path ?? '').replace(/\/+$/, '')}/v1`;
}

export interface AppContextValue {
  config: UiConfig;
  /** Shared read API (`/v4`). */
  client: ApiClient;
  /** Local-only API (`/v1`): actions, jobs, remote browse. */
  localClient: ApiClient;
}

/** Undefined only outside the provider, which is a programming error. */
export const AppContext = createContext<AppContextValue | undefined>(undefined);
