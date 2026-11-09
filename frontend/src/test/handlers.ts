import { HttpResponse, http } from 'msw';
import {
  aliasesFixture,
  assetFilesFixture,
  assetGroupsFixture,
  assetsFixture,
  emptyPage,
  genomeDetailFixture,
  genomesFixture,
  relationshipsFixture,
  remoteAssetsFixture,
  remoteGenomesFixture,
  remoteServersFixture,
  stagedAssetsFixture,
  genomeSequencesFixture,
  genomeComparisonFixture,
} from './fixtures';

/** Tests run against a relative api_base, so MSW matches path-only patterns. */
export const API = '/v4';
export const LOCAL_API = '/v1';

export const handlers = [
  http.get(`${API}/genomes`, () => HttpResponse.json(genomesFixture)),
  http.get(`${API}/genomes/:digest`, () => HttpResponse.json(genomeDetailFixture)),
  http.get(`${API}/aliases`, () => HttpResponse.json(aliasesFixture)),
  http.get(`${API}/assets`, () => HttpResponse.json(assetsFixture)),
  http.get(`${API}/assets/:digest/files`, () => HttpResponse.json(assetFilesFixture)),
  http.get(`${API}/assets/:digest`, () =>
    HttpResponse.json(assetsFixture.results[0]),
  ),
  // Asset groups list is empty by default; seed with `server.use(...)`.
  http.get(`${API}/asset_groups/:id`, () => HttpResponse.json(assetGroupsFixture.results[0])),
  http.get(`${API}/asset_groups`, () => HttpResponse.json(emptyPage)),
  http.get(`${API}/relationships/:digest`, () => HttpResponse.json(relationshipsFixture)),
  http.get(`${API}/staged_assets`, () => HttpResponse.json(stagedAssetsFixture)),
  http.get(`${API}/asset_classes`, () => HttpResponse.json(emptyPage)),
  http.get(`${API}/recipes`, () => HttpResponse.json(emptyPage)),
  http.get(`${API}/configurations`, () => HttpResponse.json(emptyPage)),
  // Server mode only (`archives`).
  http.get(`${API}/summary`, () =>
    HttpResponse.json({ genomes: 2, asset_groups: 3, assets: 4 }),
  ),
  // Jobs: empty by default.
  http.get(`${LOCAL_API}/jobs`, () => HttpResponse.json(emptyPage)),
  http.get(`${LOCAL_API}/jobs/:id/log`, () =>
    HttpResponse.json({ lines: [], next_offset: 0, has_more: false }),
  ),
  http.get(`${LOCAL_API}/remote/servers`, () => HttpResponse.json(remoteServersFixture)),
  http.get(`${LOCAL_API}/remote/genomes`, () => HttpResponse.json(remoteGenomesFixture)),
  http.get(`${LOCAL_API}/remote/assets`, () => HttpResponse.json(remoteAssetsFixture)),
  http.get(`${API}/genomes/:digest/sequences`, () => HttpResponse.json(genomeSequencesFixture)),
  http.get(`${API}/genomes/:digest/compare/:other`, () =>
    HttpResponse.json(genomeComparisonFixture),
  ),
];
