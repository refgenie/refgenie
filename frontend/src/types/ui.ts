/**
 * UI bootstrap contract: the `refgenie` object of `GET /service-info`.
 *
 * The capability key set is a shared contract across three plans (app-factory
 * emits it, this app consumes it, the manage UI reads it). Adding a key is a
 * simultaneous change in all three. A missing key reads as `false`.
 */

export type Mode = 'local' | 'server';

export interface Capabilities {
  pull: boolean;
  build: boolean;
  delete: boolean;
  aliases_write: boolean;
  subscriptions: boolean;
  recipes_write: boolean;
  asset_classes_write: boolean;
  remote_browse: boolean;
  genome_init: boolean;
  jobs: boolean;
  jobs_cancel: boolean;
  downloads: boolean;
  archives: boolean;
  seqcol: boolean;
  drs: boolean;
}

export type CapabilityKey = keyof Capabilities;

export const CAPABILITY_KEYS: readonly CapabilityKey[] = [
  'pull',
  'build',
  'delete',
  'aliases_write',
  'subscriptions',
  'recipes_write',
  'asset_classes_write',
  'remote_browse',
  'genome_init',
  'jobs',
  'jobs_cancel',
  'downloads',
  'archives',
  'seqcol',
  'drs',
] as const;

export interface UiLinks {
  docs?: string;
  github?: string;
  openapi?: string;
  [key: string]: string | undefined;
}

export interface WebUiBuildInfo {
  commit?: string | null;
  built_at?: string | null;
  dirty?: boolean | null;
}

export interface UiConfig {
  mode: Mode;
  /** Absolute or root-relative base for the shared read API, e.g. `/v4`. */
  api_base: string;
  /** Sub-path the app is mounted under; the router `basename`. */
  root_path: string;
  refgenie_version: string;
  service_name: string;
  capabilities: Capabilities;
  links: UiLinks;
  web_ui?: WebUiBuildInfo | null;
  /**
   * True when `/service-info` could not be reached and the read-only static
   * fallback is in force. The UI shows a dismissible banner in that case.
   */
  degraded: boolean;
}

/** The `refgenie` sub-object as it arrives on the wire (every field optional). */
export interface ServiceInfoRefgenieBlock {
  mode?: string;
  api_base?: string;
  root_path?: string;
  refgenie_version?: string;
  service_name?: string;
  capabilities?: Partial<Record<string, boolean>>;
  links?: UiLinks;
  web_ui?: WebUiBuildInfo | null;
}

export interface ServiceInfoResponse {
  refgenie?: ServiceInfoRefgenieBlock;
  [key: string]: unknown;
}
