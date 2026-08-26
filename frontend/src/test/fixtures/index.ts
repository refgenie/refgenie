/**
 * Typed views over the raw JSON captures.
 *
 * The JSON files are literal API responses (see frontend/README.md for the
 * capture procedure), so a shape change on the server shows up here first.
 */

import aliasesJson from './aliases.json';
import assetFilesJson from './assetFiles.json';
import assetsJson from './assets.json';
import genomeDetailJson from './genomeDetail.json';
import genomesJson from './genomes.json';
import relationshipsJson from './relationships.json';
import remoteAssetsJson from './remoteAssets.json';
import remoteGenomesJson from './remoteGenomes.json';
import remoteServersJson from './remoteServers.json';
import stagedAssetsJson from './stagedAssets.json';

import type {
  AliasPublic,
  AssetFilesResponse,
  AssetResponse,
  GenomeDetailResponse,
  GenomeResponse,
  RelationshipsExpandedResponse,
  RemoteAsset,
  RemoteGenome,
  RemoteServersResponse,
  StagedAssetPublic,
} from '../../types/api';
import type { Paginated } from '../../types/pagination';

export const genomesFixture = genomesJson as unknown as Paginated<GenomeResponse>;
export const genomeDetailFixture = genomeDetailJson as unknown as GenomeDetailResponse;
export const aliasesFixture = aliasesJson as unknown as Paginated<AliasPublic>;
export const assetsFixture = assetsJson as unknown as Paginated<AssetResponse>;
export const assetFilesFixture = assetFilesJson as unknown as AssetFilesResponse;
export const relationshipsFixture =
  relationshipsJson as unknown as RelationshipsExpandedResponse;
export const stagedAssetsFixture =
  stagedAssetsJson as unknown as Paginated<StagedAssetPublic>;
export const remoteServersFixture = remoteServersJson as unknown as RemoteServersResponse;
export const remoteGenomesFixture = remoteGenomesJson as unknown as RemoteGenome[];
export const remoteAssetsFixture = remoteAssetsJson as unknown as RemoteAsset[];

export const emptyPage = { items: [], pagination: { offset: 0, limit: 50, total: 0 } };
