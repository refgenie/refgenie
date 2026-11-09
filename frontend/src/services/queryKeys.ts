import type { JobKind } from './contracts';

/** Cache keys, each built from the same params object passed to its fetcher. */

export const qk = {
  uiConfig: () => ['ui-config'] as const,

  genomes: (p: object) => ['genomes', p] as const,
  genome: (d: string) => ['genome', d] as const,
  genomeIndex: () => ['genome-index'] as const,
  // Both start with 'genome', so the 'genome' family invalidation clears them.
  genomeSequences: (d: string) => ['genome', d, 'sequences'] as const,
  genomeComparison: (a: string, b: string) => ['genome', a, 'compare', b] as const,

  assetGroups: (p: object) => ['asset-groups', p] as const,
  assetGroup: (id: number) => ['asset-group', id] as const,

  assets: (p: object) => ['assets', p] as const,
  asset: (d: string) => ['asset', d] as const,
  assetFiles: (d: string) => ['asset-files', d] as const,
  // Local mode's own file access (`/v1/files`), not the staged-file listing.
  localAssetFiles: (d: string) => ['local-asset-files', d] as const,
  localFilePreview: (d: string, path: string) => ['local-file-preview', d, path] as const,

  assetClasses: (p: object) => ['asset-classes', p] as const,
  assetClass: (id: number) => ['asset-class', id] as const,
  assetClassIndex: () => ['asset-class-index'] as const,

  recipes: (p: object) => ['recipes', p] as const,
  recipe: (id: number) => ['recipe', id] as const,

  aliases: (p: object) => ['aliases', p] as const,

  stagedAssets: (p: object) => ['staged-assets', p] as const,
  stagedAsset: (id: number) => ['staged-asset', id] as const,

  relationships: (d: string, expand: boolean) => ['relationships', d, expand] as const,

  configurations: (p: object) => ['configurations', p] as const,
  configuration: (id: number) => ['configuration', id] as const,

  summary: () => ['summary'] as const,
  archives: (p: object) => ['archives', p] as const,

  jobs: (p: object) => ['jobs', p] as const,
  jobLog: (id: string, offset: number) => ['job-log', id, offset] as const,

  // Digests on the connected local refgenie, keyed by base URL so a new port
  // never serves a stale set.
  bridgeDigests: (baseUrl: string) => ['bridge-digests', baseUrl] as const,

  remoteServers: () => ['remote-servers'] as const,
  dataChannels: () => ['data-channels'] as const,
  remoteGenomes: (serverUrl: string | undefined) => ['remote-genomes', serverUrl] as const,
  remoteAssets: (serverUrl: string | undefined, genomeDigest: string) =>
    ['remote-assets', serverUrl, genomeDigest] as const,
};

const genomes = ['genomes', 'genome', 'genome-index', 'summary'];
const assets = [
  'assets',
  'asset',
  'asset-files',
  'local-asset-files',
  'local-file-preview',
  'asset-groups',
  'asset-group',
  'relationships',
  'staged-assets',
  'staged-asset',
];
const aliases = ['aliases', 'alias'];
const servers = ['remote-servers', 'configurations', 'configuration', 'remote-genomes', 'remote-assets'];
// A channel sync registers recipes and asset classes.
const channels = ['data-channels', 'recipes', 'recipe', 'asset-classes', 'asset-class', 'asset-class-index'];

/**
 * Key families that a finished job or a curation action makes stale. Lists are
 * refetched, never edited in place, because a refgenie change can partly fail
 * on the server.
 */
export const staleAfter = {
  pull: [...genomes, ...assets, ...aliases],
  build: [...genomes, ...assets, ...aliases],
  genome_init: [...genomes, ...aliases],
  'asset.delete': [...assets, ...genomes],
  'asset.set_default': assets,
  'genome.delete': [...genomes, ...assets, ...aliases],
  'alias.set': [...aliases, ...genomes],
  'alias.remove': [...aliases, ...genomes],
  subscribe: servers,
  unsubscribe: servers,
  'data_channel.add': channels,
  'data_channel.sync': channels,
  'data_channel.remove': ['data-channels'],
} satisfies Record<JobKind, readonly string[]> & Record<string, readonly string[]>;
