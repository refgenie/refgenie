/**
 * One function per `/v1/actions` verb (local mode only), each sent through
 * `ApiClient.mutate`. Recipe and asset-class registration stay CLI-only, since
 * they would read a server-side path or URL from the request body.
 */

import { ACTION_PATHS } from './contracts';
import type {
  ActionResult,
  AliasRequest,
  BuildRequest,
  GenomeInitRequest,
  JobRef,
  PreflightResponse,
  PullRequest,
  SetDefaultAssetRequest,
  DataChannelAddRequest,
  SubscribeRequest,
  UnsubscribeRequest,
} from './contracts';
import type { ApiClient, RequestInitLite } from './http';

// === Long operations: 202 + a JobRef ===

export const pullAsset = (c: ApiClient, body: PullRequest, init?: RequestInitLite) =>
  c.mutate<JobRef>(ACTION_PATHS.pull, {
    method: 'POST',
    action: 'pull',
    body,
    signal: init?.signal,
  });

export const buildAsset = (c: ApiClient, body: BuildRequest, init?: RequestInitLite) =>
  c.mutate<JobRef>(ACTION_PATHS.build, {
    method: 'POST',
    action: 'build',
    body,
    signal: init?.signal,
  });

export const initGenome = (c: ApiClient, body: GenomeInitRequest, init?: RequestInitLite) =>
  c.mutate<JobRef>(ACTION_PATHS.genomes, {
    method: 'POST',
    action: 'genome.init',
    body,
    signal: init?.signal,
  });

/** Check build inputs without starting a job, so the form can show errors first. */
export const preflightBuild = (c: ApiClient, body: BuildRequest, init?: RequestInitLite) =>
  c.mutate<PreflightResponse>(ACTION_PATHS.buildPreflight, {
    method: 'POST',
    action: 'build.preflight',
    body,
    signal: init?.signal,
  });

// === Short synchronous mutations: 200 + an ActionResult ===

export const deleteAsset = (c: ApiClient, digest: string, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.asset(digest), {
    method: 'DELETE',
    action: 'asset.delete',
    signal: init?.signal,
  });

/** Also deletes the genome's aliases, asset groups and assets. */
export const deleteGenome = (c: ApiClient, genomeDigest: string, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.genome(genomeDigest), {
    method: 'DELETE',
    action: 'genome.delete',
    signal: init?.signal,
  });

export const setAlias = (c: ApiClient, body: AliasRequest, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.aliases, {
    method: 'POST',
    action: 'alias.set',
    body,
    signal: init?.signal,
  });

export const removeAlias = (c: ApiClient, name: string, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.alias(name), {
    method: 'DELETE',
    action: 'alias.remove',
    signal: init?.signal,
  });

export const subscribe = (c: ApiClient, body: SubscribeRequest, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.subscriptions, {
    method: 'POST',
    action: 'subscribe',
    body,
    signal: init?.signal,
  });

export const unsubscribe = (c: ApiClient, body: UnsubscribeRequest, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.subscriptions, {
    method: 'DELETE',
    action: 'unsubscribe',
    body,
    signal: init?.signal,
  });

export const addDataChannel = (c: ApiClient, body: DataChannelAddRequest, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.dataChannels, {
    method: 'POST',
    action: 'data_channel.add',
    body,
    signal: init?.signal,
  });

export const removeDataChannel = (c: ApiClient, name: string, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.dataChannel(name), {
    method: 'DELETE',
    action: 'data_channel.remove',
    signal: init?.signal,
  });

export const syncDataChannel = (c: ApiClient, name: string, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.dataChannelSync(name), {
    method: 'POST',
    action: 'data_channel.sync',
    signal: init?.signal,
  });

export const setDefaultAsset = (
  c: ApiClient,
  body: SetDefaultAssetRequest,
  init?: RequestInitLite,
) =>
  c.mutate<ActionResult>(ACTION_PATHS.assetDefault, {
    method: 'POST',
    action: 'asset.set_default',
    body,
    signal: init?.signal,
  });
