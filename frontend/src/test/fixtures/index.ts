/** Typed views over captured API responses (see frontend/README.md). */

import aliasesJson from './aliases.json';
import assetClassesJson from './assetClasses.json';
import assetFilesJson from './assetFiles.json';
import assetGroupsJson from './assetGroups.json';
import assetsJson from './assets.json';
import genomeDetailJson from './genomeDetail.json';
import genomesJson from './genomes.json';
import recipesJson from './recipes.json';
import relationshipsJson from './relationships.json';
import remoteAssetsJson from './remoteAssets.json';
import remoteGenomesJson from './remoteGenomes.json';
import remoteServersJson from './remoteServers.json';
import stagedAssetsJson from './stagedAssets.json';
import genomeSequencesJson from './genomeSequences.json';
import genomeComparisonJson from './genomeComparison.json';
import type { GenomeSequencesResponse, SeqcolComparison } from '../../types/api';

import type {
  AliasPublic,
  AssetClassPublic,
  AssetFilesResponse,
  AssetGroupPublic,
  AssetResponse,
  GenomeDetailResponse,
  GenomeResponse,
  RecipeResponse,
  RelationshipsExpandedResponse,
  RemoteAsset,
  RemoteGenome,
  RemoteServersResponse,
  StagedAssetPublic,
} from '../../types/api';
import type { PageResponse } from '../../types/pagination';

export const genomesFixture = genomesJson as unknown as PageResponse<GenomeResponse>;
export const genomeDetailFixture = genomeDetailJson as unknown as GenomeDetailResponse;
export const aliasesFixture = aliasesJson as unknown as PageResponse<AliasPublic>;
export const assetsFixture = assetsJson as unknown as PageResponse<AssetResponse>;
export const assetGroupsFixture = assetGroupsJson as unknown as PageResponse<AssetGroupPublic>;
export const assetClassesFixture =
  assetClassesJson as unknown as PageResponse<AssetClassPublic>;
export const recipesFixture = recipesJson as unknown as PageResponse<RecipeResponse>;
export const assetFilesFixture = assetFilesJson as unknown as AssetFilesResponse;
export const relationshipsFixture =
  relationshipsJson as unknown as RelationshipsExpandedResponse;
export const stagedAssetsFixture =
  stagedAssetsJson as unknown as PageResponse<StagedAssetPublic>;
export const remoteServersFixture = remoteServersJson as unknown as RemoteServersResponse;
export const remoteGenomesFixture = remoteGenomesJson as unknown as RemoteGenome[];
export const remoteAssetsFixture = remoteAssetsJson as unknown as RemoteAsset[];

export const emptyPage = { results: [], pagination: { page: 0, page_size: 50, total: 0 } };

// === Genome sequences ===
export const genomeSequencesFixture =
  genomeSequencesJson as unknown as GenomeSequencesResponse;
export const genomeComparisonFixture =
  genomeComparisonJson as unknown as SeqcolComparison;
