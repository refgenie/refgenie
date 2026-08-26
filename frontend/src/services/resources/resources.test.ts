/**
 * One assertion per resource: the exact request path and query string produced
 * for a representative options object, and that the response comes back
 * unchanged.
 *
 * This is what proves the UI hits the server's search and pagination
 * endpoints rather than filtering a fully-loaded array on the client.
 */

import { describe, expect, it, vi } from 'vitest';
import { ApiClient } from '../http';
import { getGenome, listGenomes } from './genomes';
import { getAssetGroup, listAssetGroups } from './assetGroups';
import { getAsset, listAssetFiles, listAssets } from './assets';
import { getAssetClass, listAssetClasses } from './assetClasses';
import { getRecipe, listRecipes } from './recipes';
import { getAlias, listAliases } from './aliases';
import { listStagedAssets } from './stagedAssets';
import { getRelationships, getRelationshipsExpanded } from './relationships';
import { getConfiguration, listConfigurations } from './configurations';
import {
  archiveDownloadUrl,
  assetFileDownloadUrl,
  getSpeciesSummary,
  getSummary,
  listArchives,
} from './serverInfo';
import { listRemoteAssets, listRemoteGenomes, listRemoteServers } from './remote';

function spyClient(payload: unknown = { ok: true }) {
  const calls: string[] = [];
  const fetchImpl = vi.fn(async (url: string) => {
    calls.push(url);
    return new Response(JSON.stringify(payload), {
      headers: { 'content-type': 'application/json' },
    });
  });
  const client = new ApiClient({ baseUrl: '/v4', fetchImpl: fetchImpl as never });
  return { client, calls };
}

describe('resource request shapes', () => {
  it('genomes', async () => {
    const { client, calls } = spyClient({ items: [], pagination: {} });
    const result = await listGenomes(client, {
      q: 'hg38',
      searchFields: ['digest', 'aliases'],
      operator: 'starts_with',
      offset: 50,
      limit: 50,
    });
    expect(calls[0]).toBe(
      '/v4/genomes?q=hg38&search_fields=digest%2Caliases&operator=starts_with&offset=50&limit=50',
    );
    expect(result).toEqual({ items: [], pagination: {} });

    await getGenome(client, 'abc/def');
    expect(calls[1]).toBe('/v4/genomes/abc%2Fdef');
  });

  it('asset groups', async () => {
    const { client, calls } = spyClient();
    await listAssetGroups(client, { genome_digest: 'g1', q: 'fasta', limit: 10 });
    expect(calls[0]).toBe('/v4/asset_groups?genome_digest=g1&q=fasta&limit=10');
    await getAssetGroup(client, 3);
    expect(calls[1]).toBe('/v4/asset_groups/3');
  });

  it('assets', async () => {
    const { client, calls } = spyClient();
    await listAssets(client, { genome_digest: 'g1', q: 'fasta', searchFields: ['name'] });
    expect(calls[0]).toBe('/v4/assets?genome_digest=g1&q=fasta&search_fields=name');
    await getAsset(client, 'd1');
    expect(calls[1]).toBe('/v4/assets/d1');
    await listAssetFiles(client, 'd1');
    expect(calls[2]).toBe('/v4/assets/d1/files');
  });

  it('asset classes', async () => {
    const { client, calls } = spyClient();
    await listAssetClasses(client, { q: 'fasta', operator: 'eq' });
    expect(calls[0]).toBe('/v4/asset_classes?q=fasta&operator=eq');
    await getAssetClass(client, 9);
    expect(calls[1]).toBe('/v4/asset_classes/9');
  });

  it('recipes', async () => {
    const { client, calls } = spyClient();
    await listRecipes(client, { q: 'bwa', offset: 0 });
    expect(calls[0]).toBe('/v4/recipes?q=bwa&offset=0');
    await getRecipe(client, 7);
    expect(calls[1]).toBe('/v4/recipes/7');
  });

  it('aliases', async () => {
    const { client, calls } = spyClient();
    await listAliases(client, { q: 'hg', operator: 'contains' });
    expect(calls[0]).toBe('/v4/aliases?q=hg&operator=contains');
    await getAlias(client, 'hg 38');
    expect(calls[1]).toBe('/v4/aliases/hg%2038');
  });

  it('staged assets', async () => {
    const { client, calls } = spyClient();
    await listStagedAssets(client, { asset_digest: 'd1', mode: 'archive' });
    expect(calls[0]).toBe('/v4/staged_assets?asset_digest=d1&mode=archive');
  });

  it('relationships', async () => {
    const { client, calls } = spyClient();
    await getRelationships(client, 'd1');
    expect(calls[0]).toBe('/v4/relationships/d1?expand=false');
    await getRelationshipsExpanded(client, 'd1');
    expect(calls[1]).toBe('/v4/relationships/d1?expand=true');
  });

  it('configurations', async () => {
    const { client, calls } = spyClient();
    await listConfigurations(client, { limit: 1 });
    expect(calls[0]).toBe('/v4/configurations?limit=1');
    await getConfiguration(client, 1);
    expect(calls[1]).toBe('/v4/configurations/1');
  });

  it('server info', async () => {
    const { client, calls } = spyClient();
    await getSummary(client);
    expect(calls[0]).toBe('/v4/summary');
    await getSpeciesSummary(client);
    expect(calls[1]).toBe('/v4/species/summary');
    await listArchives(client, { genome_digest: 'g1' });
    expect(calls[2]).toBe('/v4/archives?genome_digest=g1');
  });

  it('remote browse rides the local /v1 client', async () => {
    const calls: string[] = [];
    const fetchImpl = vi.fn(async (url: string) => {
      calls.push(url);
      return new Response('[]', { headers: { 'content-type': 'application/json' } });
    });
    const local = new ApiClient({ baseUrl: '/v1', fetchImpl: fetchImpl as never });
    await listRemoteServers(local);
    expect(calls[0]).toBe('/v1/remote/servers');
    await listRemoteGenomes(local, { server_url: 'https://api.refgenie.org' });
    expect(calls[1]).toBe('/v1/remote/genomes?server_url=https%3A%2F%2Fapi.refgenie.org');
    await listRemoteAssets(local, { genome_digest: 'g1' });
    expect(calls[2]).toBe('/v1/remote/assets?genome_digest=g1');
  });
});

describe('download URL builders', () => {
  it('builds an archive download anchor target', () => {
    const client = new ApiClient({ baseUrl: 'https://api.refgenie.org/v4' });
    expect(archiveDownloadUrl(client, 'd1')).toBe(
      'https://api.refgenie.org/v4/archives/d1/download',
    );
  });

  it('encodes each path segment but keeps the / separators', () => {
    const client = new ApiClient({ baseUrl: '/v4' });
    expect(assetFileDownloadUrl(client, 'd1', 'sub dir/hg38.fa.fai')).toBe(
      '/v4/assets/d1/files/sub%20dir/hg38.fa.fai',
    );
  });
});
