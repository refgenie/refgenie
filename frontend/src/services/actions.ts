/**
 * One thin function per verb on `/v1/actions`, local mode only.
 *
 * Every one of these goes through `ApiClient.mutate`, which is the only place
 * `X-Refgenie-Action` is attached. A missing header is a 403 with code
 * `missing_action_header`, so there is no path here that can quietly skip it.
 *
 * `recipe.add` / `asset_class.add` are deliberately absent: they take a
 * server-side path or URL from an HTTP body (arbitrary local-file read plus
 * SSRF) and stay out of v1. Registration remains a CLI operation.
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

/**
 * Validates build inputs without starting a job: recipe, genome, input assets,
 * required files/params and the resolved `default_asset`. Without it the form
 * could only report failures after a job had already started and failed.
 */
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

/** Cascades unconditionally: aliases, asset groups and assets all go with it. */
export const deleteGenome = (c: ApiClient, genomeRef: string, init?: RequestInitLite) =>
  c.mutate<ActionResult>(ACTION_PATHS.genome(genomeRef), {
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
