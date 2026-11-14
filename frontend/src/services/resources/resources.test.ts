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
import { fetchAllGenomes, getGenome, listGenomes } from './genomes';
import { getAssetGroup, listAssetGroups } from './assetGroups';
import { getAsset, listAssetFiles, listAssets } from './assets';
import { getAssetClass, listAssetClasses } from './assetClasses';
import { getRecipe, listRecipes } from './recipes';
import { listAliases } from './aliases';
import { listStagedAssets } from './stagedAssets';
import { getRelationshipsExpanded } from './relationships';
import { listConfigurations } from './configurations';
import {
  archiveDownloadUrl,
  assetFileDownloadUrl,
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
    const { client, calls } = spyClient({ results: [], pagination: {} });
    const result = await listGenomes(client, {
      q: 'hg38',
      searchFields: ['digest', 'aliases'],
      operator: 'starts_with',
      page: 1,
      pageSize: 50,
    });
    expect(calls[0]).toBe(
      '/v4/genomes?q=hg38&search_fields=digest%2Caliases&operator=starts_with&page=1&page_size=50',
    );
    expect(result).toEqual({ results: [], pagination: {} });

    await getGenome(client, 'abc/def');
    expect(calls[1]).toBe('/v4/genomes/abc%2Fdef');

    await listGenomes(client, { has_assets: true, taxon_id: 9606 });
    expect(calls[2]).toBe('/v4/genomes?has_assets=true&taxon_id=9606');
  });

  it('fetches every genome page by page and stops on the last', async () => {
    const calls: string[] = [];
    const fetchImpl = vi.fn(async (url: string) => {
      calls.push(url);
      const page = Number(new URL(url, 'http://x').searchParams.get('page'));
      // 2500 rows at 1000 per page: two full pages and a short third.
      const size = page < 2 ? 1000 : 500;
      const results = Array.from({ length: size }, (_, i) => ({ digest: `g${page}-${i}` }));
      return new Response(
        JSON.stringify({ results, pagination: { page, page_size: 1000, total: 2500 } }),
        { headers: { 'content-type': 'application/json' } },
      );
    });
    const client = new ApiClient({ baseUrl: '/v4', fetchImpl: fetchImpl as never });
    const all = await fetchAllGenomes(client);
    expect(all).toHaveLength(2500);
    expect(calls).toEqual([
      '/v4/genomes?page=0&page_size=1000',
      '/v4/genomes?page=1&page_size=1000',
      '/v4/genomes?page=2&page_size=1000',
    ]);
  });

  it('asset groups', async () => {
    const { client, calls } = spyClient();
    await listAssetGroups(client, { genome_digest: 'g1', q: 'fasta', pageSize: 10 });
    expect(calls[0]).toBe('/v4/asset_groups?genome_digest=g1&q=fasta&page_size=10');
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

    // Each AssetGroup-joining filter still serializes on its own.
    await listAssets(client, { asset_group_id: 43 });
    expect(calls[3]).toBe('/v4/assets?asset_group_id=43');
    await listAssets(client, { asset_group_name: 'fasta' });
    expect(calls[4]).toBe('/v4/assets?asset_group_name=fasta');
  });

  it('rejects two AssetGroup-joining filters at compile time', async () => {
    const { client, calls } = spyClient();
    // The server joins AssetGroup once per filter, so any pair is a 500.
    // `ListAssetsParams` is a union that makes the pair unrepresentable; this
    // line fails the build if that union is ever flattened back out.
    // @ts-expect-error two AssetGroup-joining filters is a server 500; see ListAssetsParams
    await listAssets(client, { genome_digest: 'g1', asset_group_name: 'fasta' });
    expect(calls).toHaveLength(1);
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
    await listRecipes(client, { q: 'bwa', page: 0 });
    expect(calls[0]).toBe('/v4/recipes?q=bwa&page=0');
    await getRecipe(client, 7);
    expect(calls[1]).toBe('/v4/recipes/7');
  });

  it('aliases', async () => {
    const { client, calls } = spyClient();
    await listAliases(client, { q: 'hg', operator: 'contains' });
    expect(calls[0]).toBe('/v4/aliases?q=hg&operator=contains');
  });

  it('staged assets', async () => {
    const { client, calls } = spyClient();
    await listStagedAssets(client, { asset_digest: 'd1', mode: 'archive' });
    expect(calls[0]).toBe('/v4/staged_assets?asset_digest=d1&mode=archive');
  });

  it('relationships', async () => {
    const { client, calls } = spyClient();
    await getRelationshipsExpanded(client, 'd1');
    expect(calls[0]).toBe('/v4/relationships/d1?expand=true');
  });

  it('configurations', async () => {
    const { client, calls } = spyClient();
    await listConfigurations(client, { pageSize: 1 });
    expect(calls[0]).toBe('/v4/configurations?page_size=1');
  });

  it('server info', async () => {
    const { client, calls } = spyClient();
    await listArchives(client, { genome_digest: 'g1' });
    expect(calls[0]).toBe('/v4/archives?genome_digest=g1');
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
