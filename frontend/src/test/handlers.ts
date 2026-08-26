import { HttpResponse, http } from 'msw';
import {
  aliasesFixture,
  assetFilesFixture,
  assetsFixture,
  emptyPage,
  genomeDetailFixture,
  genomesFixture,
  relationshipsFixture,
  remoteAssetsFixture,
  remoteGenomesFixture,
  remoteServersFixture,
  stagedAssetsFixture,
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
    HttpResponse.json(assetsFixture.items[0]),
  ),
  http.get(`${API}/relationships/:digest`, () => HttpResponse.json(relationshipsFixture)),
  http.get(`${API}/staged_assets`, () => HttpResponse.json(stagedAssetsFixture)),
  http.get(`${API}/asset_classes`, () => HttpResponse.json(emptyPage)),
  http.get(`${API}/recipes`, () => HttpResponse.json(emptyPage)),
  http.get(`${API}/configurations`, () => HttpResponse.json(emptyPage)),
  // Jobs. Empty by default; a test that cares seeds the store or overrides
  // these with `server.use(...)`.
  http.get(`${LOCAL_API}/jobs`, () => HttpResponse.json(emptyPage)),
  http.get(`${LOCAL_API}/jobs/:id/log`, () =>
    HttpResponse.json({ lines: [], next_offset: 0, truncated: false }),
  ),
  http.get(`${LOCAL_API}/remote/servers`, () => HttpResponse.json(remoteServersFixture)),
  http.get(`${LOCAL_API}/remote/genomes`, () => HttpResponse.json(remoteGenomesFixture)),
  http.get(`${LOCAL_API}/remote/assets`, () => HttpResponse.json(remoteAssetsFixture)),
];
