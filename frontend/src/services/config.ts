/**
 * Bootstrap config resolution.
 *
 * One attempt at `/service-info`, then a read-only fallback. Never a retry
 * loop: the fallback is what keeps a pure-static, cross-origin deployment
 * possible, where `/service-info` simply does not exist on the SPA's origin.
 */

import { ApiError } from './http';
import { coerceCapabilities, noCapabilities } from './capabilities';
import type { Capabilities, Mode, ServiceInfoResponse, UiConfig } from '../types/ui';

export const SERVICE_INFO_PATH = '/service-info';

/**
 * The read-only server-mode capability set used when `/service-info` cannot be
 * reached. Reads are assumed to work; every command affordance stays off.
 */
function staticFallbackCapabilities(): Capabilities {
  return {
    ...noCapabilities(),
    downloads: true,
    archives: true,
    seqcol: true,
    drs: true,
  };
}

function coerceMode(raw: string | undefined): Mode {
  return raw === 'local' ? 'local' : 'server';
}

function envValue(name: 'VITE_SERVICE_INFO_URL' | 'VITE_API_BASE'): string | undefined {
  const value = import.meta.env?.[name];
  return typeof value === 'string' && value !== '' ? value : undefined;
}

export function fallbackConfig(): UiConfig {
  return {
    mode: 'server',
    api_base: envValue('VITE_API_BASE') ?? '/v4',
    root_path: '',
    refgenie_version: 'unknown',
    service_name: 'refgenie',
    capabilities: staticFallbackCapabilities(),
    web_ui: null,
    degraded: true,
  };
}

export function parseServiceInfo(payload: ServiceInfoResponse): UiConfig {
  const block = payload.refgenie ?? {};
  return {
    mode: coerceMode(block.mode),
    api_base: envValue('VITE_API_BASE') ?? block.api_base ?? '/v4',
    root_path: block.root_path ?? '',
    refgenie_version: block.refgenie_version ?? 'unknown',
    service_name: block.service_name ?? 'refgenie',
    capabilities: coerceCapabilities(block.capabilities),
    web_ui: block.web_ui ?? null,
    degraded: false,
  };
}

/**
 * Fetch and parse the bootstrap config. Resolves to the static fallback on any
 * failure (404, network error, non-JSON body, missing `refgenie` block).
 */
export async function loadUiConfig(fetchImpl: typeof fetch = fetch): Promise<UiConfig> {
  const url = envValue('VITE_SERVICE_INFO_URL') ?? SERVICE_INFO_PATH;
  try {
    const response = await fetchImpl(url, { headers: { accept: 'application/json' } });
    if (!response.ok) {
      throw new ApiError({
        status: response.status,
        detail: response.statusText,
        url,
      });
    }
    const payload = (await response.json()) as ServiceInfoResponse;
    if (!payload || typeof payload !== 'object' || !payload.refgenie) {
      return fallbackConfig();
    }
    return parseServiceInfo(payload);
  } catch {
    return fallbackConfig();
  }
}
