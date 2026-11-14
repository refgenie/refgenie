/**
 * Hand-written wire types for the refgenie JSON API, mirroring
 * refgenie/server/schemas.py and refgenie/db/tables.py.
 * `tests/test_ui_contract.py` checks that the keys still match.
 */

// === Enum unions (refgenie/db/tables.py) ===

export type SeekKeyType = 'file' | 'directory' | 'prefix' | 'string' | 'json';
export type StagingMode = 'file' | 'archive';

/** `InputEntities = dict[str, dict[str, Any]]` (db/tables.py). */
export type InputEntities = Record<string, Record<string, unknown>>;

// === Genomes ===

/** A genome's species. `synonyms` are used for search and the genome page only. */
export interface SpeciesResponse {
  taxon_id: number;
  scientific_name: string;
  common_name: string | null;
  synonyms: string[];
}

/** The genome list shape (no `assembly_level`). */
export interface GenomeResponse {
  digest: string;
  aliases: string[];
  description: string | null;
  asset_count: number;
  species_name: string | null;
  /** The genome's own, else its species'. */
  common_name: string | null;
  synonyms: string[] | null;
  taxon_id: number | null;
  assembly_source: string | null;
  assembly_accession: string | null;
  assembly_level: string | null;
  species: SpeciesResponse | null;
}

/** `GenomeDetailResponse` = `GenomePublic` + `taxon_uri` + `fhr` + `species`. */
export interface GenomeDetailResponse {
  digest: string;
  description: string | null;
  species_name: string | null;
  /** The genome's own, else its species'. */
  common_name: string | null;
  synonyms: string[] | null;
  taxon_id: number | null;
  assembly_source: string | null;
  assembly_accession: string | null;
  assembly_level: string | null;
  remote_url: string | null;
  taxon_uri: string | null;
  fhr: FhrMetadata | null;
  species: SpeciesResponse | null;
  /** Base URL of the remote refget store with this genome's sequences; null for local stores. */
  refget_store_url: string | null;
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

/** One seek key of an asset. */
export interface SeekKeyResponse {
  name: string;
  value: string;
  description: string | null;
  type: SeekKeyType;
  size: number | null;
}

/**
 * One name of an asset, with its build provenance. Provenance is per name because
 * two builds can share a content digest; all-null means no recorded build.
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

/** Untyped on the server: `catalog.py::list_asset_files`. */
export interface AssetFilesResponse {
  asset_digest: string;
  files: string[];
}

/** `refgenie/server/local/schemas.py::LocalFileEntry` (local mode, `/v1/files`). */
export interface LocalFileEntry {
  /** POSIX path relative to the asset directory. */
  path: string;
  type: 'file' | 'dir';
  size: number | null;
  previewable: boolean;
}

/** `GET /v1/files/assets/{asset_digest}`. */
export interface LocalAssetFiles {
  asset_digest: string;
  entries: LocalFileEntry[];
  truncated: boolean;
}

/** `GET /v1/files/assets/{asset_digest}/preview/{path}`. */
export interface LocalFilePreview {
  path: string;
  size: number;
  content: string;
  /** True when `content` is only the first `limit` bytes. */
  truncated: boolean;
  limit: number;
}

/** Untyped on the server: `catalog.py::get_asset_relationships` (expand=false). */
export interface RelationshipsResponse {
  asset_digest: string;
  parents: string[];
  children: string[];
}

/** Untyped on the server: `catalog.py::get_asset_relationships` (expand=true). */
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

/** One seek key an asset class declares (`AssetClassSeekKeyResponse`). */
export interface AssetClassSeekKeyResponse {
  name: string;
  value: string | null;
  description: string | null;
  type: SeekKeyType;
}

/** `GET /asset_classes/{id}`: one asset class with its seek keys. */
export interface AssetClassResponse extends AssetClassPublic {
  seek_keys: AssetClassSeekKeyResponse[];
}

export interface RecipeResponse {
  id: number | null;
  name: string;
  version: string;
  description: string | null;
  output_asset_class_id: number;
  output_asset_class_name: string | null;
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

/** One data channel the local instance syncs recipes and asset classes from. */
export interface LocalDataChannel {
  name: string;
  type: string;
  index_address: string;
  description: string | null;
  credentials_set: boolean;
  /** Always false today; the field is where trust will be read from. */
  trusted: boolean;
}

export interface LocalDataChannelsResponse {
  channels: LocalDataChannel[];
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

// === Genome sequences ===
// `GET /v4/genomes/{digest}/sequences` and `/compare/{other}` (`catalog.py`).

export interface AliasEntry {
  namespace: string;
  alias: string;
}

export interface SequenceResponse {
  name: string;
  length: number;
  /** Bare GA4GH digest, no `SQ.` prefix. */
  sha512t24u: string;
  md5: string | null;
  alphabet: string | null;
  description: string | null;
  aliases: AliasEntry[];
}

export interface GenomeSequencesResponse {
  digest: string;
  n_sequences: number;
  total_length: number;
  /** Level-1 attribute digests: names, lengths, sequences, ... */
  level1: Record<string, string>;
  collection_aliases: AliasEntry[];
  sequences: SequenceResponse[];
}

/** The GA4GH seqcol comparison (`refget.utils.compare_seqcols`). */
export interface SeqcolComparison {
  digests: { a: string; b: string };
  attributes: { a_only: string[]; b_only: string[]; a_and_b: string[] };
  array_elements: {
    a_count: Record<string, number>;
    b_count: Record<string, number>;
    a_and_b_count: Record<string, number>;
    a_and_b_same_order: Record<string, boolean | null>;
  };
  /**
   * refgenie's `/genomes/{a}/compare/{b}` only: one sentence on how the two
   * differ, ready to show. Absent from older servers.
   */
  summary?: string;
  /** With `summary`: A's sequence name -> B's, for sequences B has under another name. */
  name_map?: Record<string, string>;
}
