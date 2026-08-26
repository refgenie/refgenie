# refgenie web UI

A single React + TypeScript SPA over the refgenie JSON API. It serves **both**
the local dash (`refgenie dash`) and the public server (`refgenie serve`);
server mode is local mode with the command surface switched off.

The UI never branches on the mode string. Every affordance is gated on a
**capability flag** from `GET /service-info`, so a future authenticated server
can flip individual flags without a frontend change.

## Layout

```
src/
  app/          router table, nav registry, config provider, query client
  components/   common/ layout/ genomes/ assets/ assetClasses/ recipes/ remote/
  pages/        one file per route
  services/     ApiClient, config loader, query keys, resources/ (one per endpoint family)
  hooks/        useApiClient, useUiConfig, useCapability, useSearchParamsState, queries/
  types/        api.ts (wire types), ui.ts (bootstrap contract), pagination.ts
  styles/       tokens.css, utilities.css, components.css, modal.css, main.css
  test/         MSW setup, fixtures, renderWithProviders
```

## Development

```bash
# terminal 1 (repo root)
uv run refgenie dash            # local mode, port 8080
#   or: uv run refgenie serve   # server mode, port 8000

# terminal 2
cd frontend && npm install
npm run dev                     # http://localhost:5173
```

The Vite dev server proxies `/v4`, `/v1`, `/service-info`, `/ping`,
`/openapi.json`, `/seqcol`, `/ga4gh` and `/data_channel` to
`http://127.0.0.1:8080`. Point it elsewhere with
`VITE_DEV_BACKEND=http://127.0.0.1:8000` (see `.env.development`). Dev and
production therefore use identical relative URLs.

To develop against a cross-origin API with no `/service-info` on this origin,
set `VITE_API_BASE`; the UI falls back to a read-only server-mode config and
shows a banner.

## Scripts

| Script | What it does |
|---|---|
| `npm run dev` | Vite dev server on 5173 |
| `npm run build` | `tsc -b` then `vite build` into `../refgenie/server/webui` |
| `npm run typecheck` | TypeScript only |
| `npm run lint` | ESLint |
| `npm run test` | Vitest (jsdom + MSW) |
| `npm run style-guard` | Style-system enforcement (see below) |

## Build output

`npm run build` writes to **`../refgenie/server/webui/`**, with the hashed
bundle under `webui/_app/`. `_app` is load-bearing: Vite's default `assets`
would collide with the SPA's own `/assets` browse route and with the
`/v4/assets` API family.

`index.html` ships a literal `<base href="/">`. The server rewrites that one
token at app construction for sub-path deployments — do not remove it and do
not template it. The router `basename` and the API base come from
`/service-info` at runtime, so there is no per-deployment rebuild.

The bundle is gitignored and never committed.

## Styling rules

The style system comes from the org Web Design Style guide.

- **`src/styles/utilities.css` and `src/styles/modal.css` are verbatim copies
  from the skill. Never edit them.** `npm run style-guard` compares their
  SHA-256 against a recorded hash. If the skill itself changes, re-copy the
  file and run `npm run style-guard -- --update`.
- **`src/styles/tokens.css` is the only file allowed to contain raw values.**
  Every hex color, px and rem literal lives there; everything else references
  a `var(--…)`.
- **`src/styles/components.css` is BEM only** (`rg-block__element--modifier`).
  Reach for a utility class first; add a block only when a utility composition
  cannot express the rule.
- **No inline styles, no Tailwind, no Bootstrap, no CSS-in-JS, no web-font
  `@import`.** The local dash must render fully offline. ESLint fails on a JSX
  `style` attribute; `style-guard` fails on framework imports and `!important`.
- All modals use `src/components/common/BaseModal.tsx` (also a verbatim copy).
  No inline modal JSX.

## Testing

Vitest + jsdom + Testing Library, with MSW v2 as the network layer. An
unhandled request fails the test: if the UI hits an endpoint nobody declared,
that is a bug.

```bash
npm run test
npm run test:watch
```

### Regenerating the MSW fixtures

`src/test/fixtures/*.json` are literal API responses. Recapture them from a
real instance so the tests track reality:

```bash
uv run refgenie init            # once
uv run refgenie pull <genome>/<asset>
uv run refgenie dash &          # port 8080

curl -s localhost:8080/v4/genomes                 > src/test/fixtures/genomes.json
curl -s localhost:8080/v4/genomes/<digest>        > src/test/fixtures/genomeDetail.json
curl -s localhost:8080/v4/aliases                 > src/test/fixtures/aliases.json
curl -s "localhost:8080/v4/assets?genome_digest=<digest>" > src/test/fixtures/assets.json
curl -s localhost:8080/v4/assets/<asset>/files    > src/test/fixtures/assetFiles.json
curl -s "localhost:8080/v4/relationships/<asset>?expand=true" > src/test/fixtures/relationships.json
curl -s "localhost:8080/v4/staged_assets?asset_digest=<asset>" > src/test/fixtures/stagedAssets.json
curl -s localhost:8080/v1/remote/servers          > src/test/fixtures/remoteServers.json
curl -s localhost:8080/v1/remote/genomes          > src/test/fixtures/remoteGenomes.json
curl -s "localhost:8080/v1/remote/assets?genome_digest=<digest>" > src/test/fixtures/remoteAssets.json
```

`src/test/fixtures/index.ts` types the captures against `src/types/api.ts`.

## Keeping the wire types honest

`src/types/api.ts` is hand-written, not generated. The Python test
`tests/test_web.py` reads the app's OpenAPI schema and asserts that each
response model still has exactly the properties this file declares, so a
backend rename fails a fast unit test that names the field to change.

## The management surface

Everything that *executes* a refgenie operation lives behind a capability flag
and reports into one place.

```
src/
  services/contracts.ts   the /v1 endpoint paths and job/action wire types
  services/actions.ts     one function per action verb
  services/jobs.ts        job reads + the single SSE URL
  stores/jobStore.ts      live jobs, per-job log ring buffer, connection state
  hooks/useJobEvents.ts   the one event transport (SSE, with a polling fallback)
  components/jobs/        console, card, progress, log tail, error block, detail modal
  components/actions/     pull, delete asset, delete genome, set default, init genome
  components/build/       the recipe-driven build form
  components/manage/      subscriptions, aliases, recipes, asset classes
  pages/                  BuildPage, ManagePage, JobsPage
```

Rules that are easy to break and expensive to debug:

- **Mutations go through `ApiClient.mutate`.** It is the only place the
  `X-Refgenie-Action` header is attached, and a missing header is a 403. An
  ESLint rule stops a component calling `fetch` directly.
- **One SSE stream, never one per job.** Browsers cap ~6 HTTP/1.1 connections
  per origin. `useJobEvents` is mounted once, by `AppLayout`.
- **`EventSource` cannot send headers**, so `/v1/jobs/events` is a plain GET and
  must stay that way.
- **A duplicate submission is a 202 with `duplicate: true`, not a 409.** The UI
  focuses the existing job card; there is no error branch to write.
- **`force` is always sent explicitly on a pull.** Omitting it lets the server
  reach a stdin prompt that hangs a worker forever.
- **Refetch after an action; never edit a list optimistically.** refgenie
  mutations touch the filesystem *and* the database and have rollback paths. The
  invalidation bus (`services/invalidation.ts`) publishes domain keys, and
  `hooks/useInvalidationBridge.ts` maps them to query keys.

Recipe and asset-class *registration* are deliberately absent: they take a
server-side path or URL from an HTTP body, so they stay CLI-only in v1 and both
panels are read-only.

## Manual smoke test

The automated tests never touch a backend, so the transport gets one manual
pass:

```bash
# terminal 1 (repo root)
uv run refgenie dash            # local mode, port 8080

# terminal 2
cd frontend && npm run dev      # http://localhost:5173
```

1. Subscribe to `https://refgenomes.databio.org` on `/manage`.
2. Open `/remote`, pick a genome, and pull an asset. The job console should
   appear at the bottom of the window and update live; the header indicator
   should read **live**.
3. Kill the backend mid-pull. The indicator should go **offline**, then
   **polling**, and recover to **live** when you start the backend again — with
   no progress lost, because the reconnect resumes from the event cursor.
4. Start a build from `/build` and confirm the log tail fills in. Builds have no
   percentage, so the tail is the whole signal.
5. Restart the backend and open `/jobs`. An empty history is correct: jobs live
   in the server process and are not persisted.
