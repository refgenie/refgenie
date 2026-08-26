/**
 * Wire types for the refgenie JSON API.
 *
 * These are hand-written and mirror the pydantic / SQLModel definitions
 * field-for-field. `tests/test_ui_contract.py` in the Python package asserts
 * that the OpenAPI component schemas still have exactly these keys, so a
 * backend rename fails a fast unit test naming this file.
 *
 * Sources: refgenie/server/schemas.py and refgenie/db/tables.py.
 */

// === Enum unions (refgenie/db/tables.py) ===

export type ServingMode = 'file' | 'archive' | 'none';
export type SeekKeyType = 'file' | 'directory' | 'prefix' | 'string' | 'json';
export type RemoteType = 's3' | 'http' | 'https';
export type StagingMode = 'file' | 'archive';

/** `InputEntities = dict[str, dict[str, Any]]` (db/tables.py). */
export type InputEntities = Record<string, Record<string, unknown>>;

// === Genomes ===

/** `GenomeResponse` — the list shape. `assembly_level` is deliberately absent. */
export interface GenomeResponse {
  digest: string;
  aliases: string[];
  description: string | null;
  asset_count: number;
  species_name: string | null;
  common_name: string | null;
  taxon_id: number | null;
  assembly_source: string | null;
  assembly_accession: string | null;
}

/** `GenomeDetailResponse` = `GenomePublic` + `taxon_uri` + `fhr`. */
export interface GenomeDetailResponse {
  digest: string;
  description: string | null;
  species_name: string | null;
  common_name: string | null;
  taxon_id: number | null;
  assembly_source: string | null;
  assembly_accession: string | null;
  assembly_level: string | null;
  remote_url: string | null;
  taxon_uri: string | null;
  fhr: FhrMetadata | null;
}

// === Asset groups ===

export interface AssetGroupPublic {
  id: number | null;
  name: string;
  description: string | null;
  genome_digest: string;
  asset_class_id: number;
}

// === Assets ===

/** `SeekKeyResponse` — one seek key of an asset. */
export interface SeekKeyResponse {
  name: string;
  value: string;
  description: string | null;
  type: SeekKeyType;
  size: number | null;
}

/**
 * `AssetNameResponse` — one name of an asset, with the build behind it.
 * Provenance is per-name: two builds can share one content digest, and each
 * records its own. All-null provenance means the name records no build.
 */
export interface AssetNameResponse {
  name: string;
  is_default: boolean;
  build_digest: string | null;
  build_level1: Record<string, string> | null;
  build_digest_scheme: string | null;
  build_timestamp: string | null;
  refgenie_version: string | null;
  inputs: Record<string, unknown> | null;
  docker_image: string | null;
  docker_image_digest: string | null;
  recipe_id: number | null;
}

/** `AssetResponse` = `AssetPublic` + resolved extras. */
export interface AssetResponse {
  digest: string | null;
  name: string;
  description: string | null;
  recipe_id: number | null;
  asset_group_id: number | null;
  size: number | null;
  serving_modes_override: string[] | null;
  colocate: Array<Record<string, string>> | null;
  serving_modes: string[] | null;
  asset_class_name: string | null;
  asset_group_name: string | null;
  genome_digest: string | null;
  names: AssetNameResponse[] | null;
  seek_keys: SeekKeyResponse[] | null;
  /** Whether this asset is its group's default (any of its names is flagged). */
  is_default?: boolean | null;
}

/** Untyped on the server: `shared.py::list_asset_files`. */
export interface AssetFilesResponse {
  asset_digest: string;
  files: string[];
}

/** Untyped on the server: `shared.py::get_asset_relationships` (expand=false). */
export interface RelationshipsResponse {
  asset_digest: string;
  parents: string[];
  children: string[];
}

/** Untyped on the server: `shared.py::get_asset_relationships` (expand=true). */
export interface RelationshipsExpandedResponse {
  asset_digest: string;
  parents: AssetResponse[];
  children: AssetResponse[];
}

// === Asset classes and recipes ===

export interface AssetClassPublic {
  id: number | null;
  name: string;
  version: string;
  description: string | null;
  serving_modes: string[];
}

export interface RecipePublic {
  id: number | null;
  name: string;
  version: string;
  description: string | null;
  output_asset_class_id: number;
  command_templates: string[];
  input_params: InputEntities | null;
  input_files: InputEntities | null;
  input_assets: InputEntities | null;
  docker_image: string | null;
  custom_seek_keys: Record<string, string> | null;
  default_asset: string;
  inherent: string[] | null;
}

// === Staging, aliases, configuration ===

export interface StagedAssetPublic {
  asset_digest: string;
  mode: StagingMode;
  directory_contents: string[];
  build_commands: string[];
  download_count: number;
  tarball_digest: string | null;
  tarball_size: number | null;
}

export interface AliasPublic {
  name: string;
  genome_digest: string;
}

export interface AliasResponse {
  alias: string;
  digest: string;
  source: string;
  collection: Record<string, unknown>;
  fhr: FhrMetadata | null;
}

export interface ConfigurationPublic {
  version: number;
  servers: string[];
  genome_folder: string;
  genome_stage_folder: string | null;
}

// === Server-mode only (refgenie/server/routers/version4.py) ===

export interface ArchiveRecord {
  digest: string;
  asset_digest: string;
  tarball_digest: string | null;
  size: number | null;
  directory_contents: string[] | null;
  build_commands: string[] | null;
  download_count: number;
}

export interface DatabaseSummaryResponse {
  genomes: number;
  asset_groups: number;
  assets: number;
}

export interface SpeciesStatistics {
  genomes: number;
  asset_classes: number;
  assets: number;
}

export type SpeciesSummaryResponse = Record<string, SpeciesStatistics>;

// === FHR sidecar ===

export interface FhrPerson {
  name?: string | null;
  [key: string]: unknown;
}

export interface FhrVitalStats {
  N50?: number | null;
  L50?: number | null;
  L90?: number | null;
  totalBasePairs?: number | null;
  numberContigs?: number | null;
  numberScaffolds?: number | null;
  readTechnology?: string | null;
  [key: string]: unknown;
}

/**
 * The FHR sidecar blob. Only the fields the UI renders are named; the rest is
 * carried by the index signature and shown in the raw-JSON toggle.
 */
export interface FhrMetadata {
  vitalStats?: FhrVitalStats | null;
  metadataAuthor?: FhrPerson[] | null;
  assemblyAuthor?: FhrPerson[] | null;
  relatedLink?: string[] | null;
  dateCreated?: string | null;
  license?: string | null;
  scholarlyArticle?: string | null;
  funding?: string | null;
  accessionID?: { url?: string | null; [key: string]: unknown } | null;
  version?: string | null;
  [key: string]: unknown;
}

// === Remote browse (local mode only, /v1/remote/*) ===

export interface RemoteServer {
  url: string;
  subscribed: boolean;
  reachable: boolean;
  error: string | null;
}

export interface RemoteServersResponse {
  servers: RemoteServer[];
}

export interface RemoteGenome {
  server_url: string;
  genome_digest: string;
  aliases: string[];
  description: string | null;
}

export interface RemoteAsset {
  server_url: string;
  genome_digest: string;
  asset_group_name: string;
  asset_name: string;
  archive_size: number | null;
  asset_digest: string;
  archive_digest: string | null;
}
