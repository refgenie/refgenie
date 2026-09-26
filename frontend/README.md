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
  app/          router table, nav registry, config provider, resource cache policy
  components/   common/ layout/ genomes/ assets/ assetClasses/ recipes/ remote/ bridge/
  pages/        one file per route
  services/     ApiClient, config loader, resource cache, cache keys, resources/ (one per endpoint family), bridge/
  hooks/        useResource, useApiClient, useUiConfig, useCapability, useSearchParamsState, queries/
  types/        api.ts (wire types), ui.ts (bootstrap contract), pagination.ts
  styles/       tokens.css, utilities.css, components.css, modal.css, main.css
  test/         MSW setup, fixtures, renderWithProviders
```

### The page head

Every inner page opens with `components/layout/MiniHero` — the landing hero at
inner-page scale, carrying the breadcrumb, the `<h1>`, an optional explainer and
an optional actions cluster. **It owns the document title**, so a page must
render it in every state, passing `documentTitle={false}` while its record is
still loading; `NotAvailablePage` sets its own title because gated pages
early-return it instead. The explainer rule: a **list** page defines the concept
it lists, once, in plain words; a **detail** page carries one only where no list
page above it defines the word, which today is `/assets/:digest` and
`/asset-groups/:id` and nothing else. `LandingPage` keeps the full `.rg-hero`
and is not a MiniHero page.

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
would collide with the SPA's own `/assets/:digest` route and with the
`/v4/assets` API family.

`index.html` ships a literal `<base href="/">`. The server rewrites that one
token at app construction for sub-path deployments — do not remove it and do
not template it. The router `basename` and the API base come from
`/service-info` at runtime, so there is no per-deployment rebuild.

The bundle is gitignored and never committed.

Everything in `public/` ships at the bundle root, including `SKILL.md` — the
agent-facing capability doc served at `/SKILL.md`. Note that a *missing* file
there does not 404: the SPA catch-all answers with the shell, so Vite, the
Python server and Cloudflare all return `200 text/html`. Assert on
`Content-Type`, never on status.

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
`tests/server/test_app.py` reads the app's OpenAPI schema and asserts that each
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
  `hooks/useInvalidationBridge.ts` maps them to cache-key prefixes.

Recipe and asset-class *registration* are deliberately absent: they take a
server-side path or URL from an HTTP body, so they stay CLI-only in v1 and both
panels are read-only.

## The localhost bridge

The public site can talk to a `refgenie dash` running on the visitor's own
machine: probe it, badge the assets they already have, browse it under
`/local/*`, and pull to it.

```
src/
  services/bridge/        the /ping contract, the probe, persistence, clients, pull, digests
  stores/bridgeStore.ts   the seven connection statuses
  hooks/useBridge.ts      store read + connect/disconnect (no effect)
  hooks/useBridgeAutoConnect.ts   the ONE remembered re-probe, mounted by AppLayout
  components/bridge/      control, connect card, connect dialog, /local scope, pull button
```

Properties that are load-bearing, not incidental:

- **The probe touches exactly one port and one endpoint, never a scan.** A page
  enumerating ports on a visitor's machine is the attack this must stay
  distinguishable from. `services/bridge/probe.ts` is also the one place in the
  app that calls `fetch` outside `ApiClient`, because its correctness *is* its
  request shape: a CORS simple request, no custom headers, `credentials: 'omit'`.
- **The first probe is user-initiated.** Only a remembered `{port, instanceId}`
  licenses a later silent re-probe, which is what puts Chrome's Local Network
  Access prompt in a context the user understands.
- **`/ping` is validated, never trusted.** Any program can listen on that port,
  so validation is a collision check and the mitigation is presentational:
  local data always lives behind `/local/*` with a banner, and never merges into
  a remote list without a badge.
- **`/local/*` is read-only by construction.** `LocalScope` masks every command
  capability to false regardless of what the ping claimed; the only cross-origin
  mutation that exists anywhere is `POST /v1/actions/pull`.
- **Pull availability comes from the ping, not the page.** On refgenie.org the
  page's own `pull` capability is false, so `useCapability('pull')` would
  silently kill the feature.
- **WebKit is not worked around.** It blocks HTTPS→localhost outright (bug
  171934), so the doomed request is never issued and the card explains why.
- **The two `localStorage` keys are frozen.** `refgenie.localBridge` and
  `refgenie.localBridge.safariDismissed` carry a returning visitor's connection
  across an SPA swap on the same origin.

Adding a `/local` route means editing `SPA_CLIENT_ROUTES` in **both**
`src/app/routes.tsx` and `refgenie/server/const.py`;
`tests/server/test_app.py::test_spa_route_mirror` fails on any drift.

### Testing the bridge end to end

The bridge needs **two origins**: a page in server mode, and a `refgenie dash`
it can reach. The default dev loop above gives you neither, so this is its own
setup.

Two things trip people up before anything else does:

- **The bridge UI only mounts when `/service-info` reports `mode !== 'local'`**
  (`AppLayout.tsx`). Point `npm run dev` at a `refgenie dash` — the documented
  default — and every bridge control is correctly hidden, because the page
  thinks it *is* the local refgenie.
- **`refgenie/server/webui/` must exist.** A source checkout has no bundle and
  `refgenie dash` answers `/` with a 503. Run `npm run build` after every
  frontend change; `build.emptyOutDir` replaces the directory wholesale.

```bash
# terminal 1 — the LOCAL side (what the bridge connects TO)
env -u REFGENIE \
  REFGENIE_BRIDGE_ORIGINS='https://refgenie.org,http://localhost:5173,http://127.0.0.1:5173' \
  uv run refgenie dash --bridge full

# terminal 2 — the PUBLIC page, in server mode, over the real public API
VITE_SERVICE_INFO_URL=https://api.refgenie.org/service-info \
VITE_API_BASE=https://api.refgenie.org/v4 \
npm run dev
```

(`env -u REFGENIE` is only for machines migrating from refgenie 0.x, where a
leftover `$REFGENIE` export makes every command print a legacy-migration
warning. Drop it if you never set that variable.)

**`refgenie dash` dies with `TypeError: CORSMiddleware.__init__() got an
unexpected keyword argument 'allow_private_network'`?** Your virtualenv predates
`uv.lock`. That keyword needs Starlette 1.4.1 / FastAPI 0.141.1, and
`install_local_security()` passes it on every `dash` run regardless of
`--bridge`, so local mode cannot start at all on an older pair. Run `uv sync`.
Server mode (`refgenie serve`) is unaffected — it installs no CORS middleware.

Open `http://localhost:5173` in **Chrome**. Safari cannot test this: WebKit
blocks HTTPS→localhost outright (bug 171934) and `isWebKitOnly()` skips the
doomed request on purpose.

Both `VITE_*` variables are load-bearing. Drop `VITE_SERVICE_INFO_URL` and the
page reads `/service-info` from the dash through the proxy and hides the bridge.
Drop `VITE_API_BASE` and the remote document's relative `api_base: "/v4"`
resolves against `localhost:5173`, so the proxy points the "remote" pane at your
own dash.

Two flags on the dash side, and they do different jobs:

| Flag / variable | What it unlocks |
| --- | --- |
| `--bridge read` (default) | Cross-origin `/ping`, `/v4`, `/v1/jobs`. Enough for connect, badges and `/local/*` browsing. |
| `--bridge full` | Additionally `POST /v1/actions/pull`. **"Pull to my refgenie" is disabled without it** — `BridgePullButton` gates on `ping.bridge.actions_cross_origin`, so `read` mode shows the disabled button plus the deep-link fallback. |
| `REFGENIE_BRIDGE_ORIGINS` | The exact origins admitted. Comma-separated, and it **replaces** the default rather than appending — repeat the production origin. |

`http://localhost:5173` is already in `LocalSecuritySettings.allowed_origins`,
so the connect step works without any env var. Listing it in
`REFGENIE_BRIDGE_ORIGINS` anyway is what makes the test faithful:
`require_action_origin()` classifies bridge origins *before* trusted dev
origins, so this routes the pull through the same `bridge_mode` gate and
`BRIDGE_CROSS_ORIGIN_ACTIONS` path allowlist that `refgenie.org` will hit.

**Connecting reports "absent"?** Check which port Vite actually bound. There is
no `strictPort`, so a busy 5173 silently becomes 5174 — which is on no
allowlist, and the browser collapses that CORS rejection, a Local Network Access
denial and connection-refused into one indistinguishable `TypeError`. Free 5173
or add the port you got.

**`REFGENIE_LOCAL_ALLOWED_ORIGINS` is not the variable you want.**
`allowed_origins` is a plain `list[str]`, so it takes JSON
(`'["http://localhost:5174"]'`) and replaces the defaults, while
`bridge_origins` is `NoDecode` with a comma-splitting validator. Use
`REFGENIE_BRIDGE_ORIGINS`.

For dev servers and Cloudflare preview deployments on unpredictable origins,
`REFGENIE_BRIDGE_ORIGIN_REGEX` exists — e.g. `http://localhost:51[0-9]{2}`.
`install_local_security()` refuses `.*`-class patterns at construction. Prefer
exact origins.

#### Without a browser

Every server-side control answers to curl, which is the fastest way to tell a
policy rejection from a page bug:

```bash
curl -si -H 'Origin: https://refgenie.org' http://localhost:8080/ping
#   200, access-control-allow-origin echoed, cache-control: no-store

curl -si -H 'Origin: https://evil.example' http://localhost:8080/ping \
  | grep -i access-control-allow-origin
#   no output — no CORS grant

curl -si -H 'Host: evil.example' http://127.0.0.1:8080/ping | head -3
#   421 forbidden_host (the DNS-rebinding guard)

curl -si -X OPTIONS http://localhost:8080/v1/actions/pull \
  -H 'Origin: https://refgenie.org' \
  -H 'Access-Control-Request-Method: POST' \
  -H 'Access-Control-Request-Headers: content-type,x-refgenie-action' \
  -H 'Access-Control-Request-Private-Network: true'
#   200 with access-control-allow-private-network: true

curl -si -X POST http://localhost:8080/v1/actions/pull \
  -H 'Content-Type: application/json' \
  -d '{"asset_group":"fasta","genome":"rCRSd","force":false}' | head -10
#   403 missing_action_header (the CSRF control)
```

#### Against the deployed site

```bash
env -u REFGENIE \
  REFGENIE_BRIDGE_ORIGINS='https://refgenie.org,https://api.refgenie.org' \
  uv run refgenie dash --bridge full
```

The SPA currently deploys to `https://api.refgenie.org`, which is **not** in the
shipped allowlist — hence the extra origin. Check
`https://api.refgenie.org/service-info` → `refgenie.web_ui.commit` before
concluding a missing control is a bug; the deployed bundle may predate your
branch.

#### Confirming the dash serves *your* bundle

`refgenie dash` deliberately shows no bridge UI — that page *is* the local
refgenie (`mode === 'local'`, plus `isSelfLocal()` as a second gate). Use it to
check the bundle instead:

```bash
uv run refgenie dash
curl -s http://localhost:8080/service-info | python3 -m json.tool | grep -A4 web_ui
```

`built_at` should be minutes old and `commit` should match your `git rev-parse
--short HEAD`. A 503 on `/` means `npm run build` has not run.

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
