/**
 * UI bootstrap config: the `refgenie` object of `GET /service-info`. Capability
 * keys must match what the server emits; a missing key reads as `false`.
 */

export type Mode = 'local' | 'server';

export interface Capabilities {
  pull: boolean;
  build: boolean;
  delete: boolean;
  aliases_write: boolean;
  subscriptions: boolean;
  data_channels: boolean;
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
  /** View and download asset files from disk (`/v1/files`). Never true via the bridge. */
  local_files: boolean;
}

export type CapabilityKey = keyof Capabilities;

export const CAPABILITY_KEYS: readonly CapabilityKey[] = [
  'pull',
  'build',
  'delete',
  'aliases_write',
  'subscriptions',
  'data_channels',
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
  'local_files',
] as const;

export interface WebUiBuildInfo {
  commit?: string | null;
  built_at?: string | null;
  dirty?: boolean | null;
}

export interface CollectionLink {
  label: string;
  url: string;
}

/** A per-genome link; placeholders `{assembly_accession}`, `{species_name}`, `{species_underscored}`, `{digest}`. */
export interface GenomeLinkTemplate {
  label: string;
  template: string;
}

export type CollectionColumn =
  | 'assembly_accession'
  | 'haplotype'
  | 'assembly_level'
  | 'assembly_source';

/** A branded slice of the catalog drawn from one or more stores (server config). */
export interface CollectionConfig {
  id: string;
  name: string;
  short_name: string | null;
  blurb: string;
  note: string | null;
  accent_color: string;
  logo_url: string | null;
  stores: string[];
  links: CollectionLink[];
  genome_links: GenomeLinkTemplate[];
  columns: CollectionColumn[];
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
  web_ui?: WebUiBuildInfo | null;
  /** Collection views the server declares; empty in local mode. */
  collections: CollectionConfig[];
  /** True when `/service-info` was unreachable and the read-only fallback is in use. */
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
  // The `links` block is ignored; project links live in services/projectLinks.ts.
  web_ui?: WebUiBuildInfo | null;
  collections?: unknown;
}

export interface ServiceInfoResponse {
  refgenie?: ServiceInfoRefgenieBlock;
  [key: string]: unknown;
}
