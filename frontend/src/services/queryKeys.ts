/**
 * Query keys. Every key is derived from the same options object passed to the
 * fetcher, so keys are stable and reviewable next to the request they cache.
 */

export const qk = {
  uiConfig: () => ['ui-config'] as const,

  genomes: (p: object) => ['genomes', p] as const,
  genome: (d: string) => ['genome', d] as const,

  assetGroups: (p: object) => ['asset-groups', p] as const,
  assetGroup: (id: number) => ['asset-group', id] as const,

  assets: (p: object) => ['assets', p] as const,
  asset: (d: string) => ['asset', d] as const,
  assetFiles: (d: string) => ['asset-files', d] as const,

  assetClasses: (p: object) => ['asset-classes', p] as const,
  assetClass: (id: number) => ['asset-class', id] as const,

  recipes: (p: object) => ['recipes', p] as const,
  recipe: (id: number) => ['recipe', id] as const,

  aliases: (p: object) => ['aliases', p] as const,
  alias: (name: string) => ['alias', name] as const,

  stagedAssets: (p: object) => ['staged-assets', p] as const,
  stagedAsset: (id: number) => ['staged-asset', id] as const,

  relationships: (d: string, expand: boolean) => ['relationships', d, expand] as const,

  configurations: (p: object) => ['configurations', p] as const,
  configuration: (id: number) => ['configuration', id] as const,

  summary: () => ['summary'] as const,
  speciesSummary: () => ['species-summary'] as const,
  archives: (p: object) => ['archives', p] as const,

  // Jobs live in the server process and are not persisted; these keys exist so
  // the history view participates in the same invalidation bus as everything
  // else, not because the records survive a restart.
  jobs: (p: object) => ['jobs', p] as const,
  job: (id: string) => ['job', id] as const,
  jobLog: (id: string, offset: number) => ['job-log', id, offset] as const,

  remoteServers: () => ['remote-servers'] as const,
  remoteGenomes: (serverUrl: string | undefined) => ['remote-genomes', serverUrl] as const,
  remoteAssets: (serverUrl: string | undefined, genomeDigest: string) =>
    ['remote-assets', serverUrl, genomeDigest] as const,
};
