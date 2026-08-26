# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Federated multi-store serving. A single refgenie server can serve genomes
  from several physically separate refget stores at once. A new `store`
  registry table is the single source of truth for what the service serves,
  managed with `refgenie store add/sync/list/remove/conflicts`. Content is
  digest-addressed, so identical genomes across stores dedup automatically; the
  only true conflict -- two stores mapping the same alias to different
  collection digests -- is resolved by store priority (lower integer wins), and
  the losing mapping stays reachable as a qualified `store::alias` name.
  Genomes carry a `store_name` owner; a `RefgetStoreRouter` dispatches reads to
  the owning store. Servers register their stores on boot from the
  `REFGENIE_STORES` environment variable (see `deployment/stores.example.json`).
- Background jobs in local mode. `refgenie dash` runs pull and build on
  dedicated worker threads instead of blocking a request handler, and reports
  status, progress, log lines, result and error per job:
  `POST /v1/jobs`, `GET /v1/jobs`, `GET /v1/jobs/{id}`, `GET /v1/jobs/{id}/log`,
  `POST /v1/jobs/{id}/cancel`, `DELETE /v1/jobs/{id}`, one multiplexed
  server-sent-events stream at `GET /v1/jobs/events`, and a polling fallback at
  `GET /v1/jobs/{id}/events?format=json`. Jobs live in one process's memory and
  do not survive a restart, which is why local mode runs a single worker.
- `refgenie.progress`: an optional, stdlib-only progress sink for long-running
  library operations. With no sink installed -- the CLI's situation always --
  `emit()` is a no-op and nothing changes. Downloads report byte counts to it
  when one is installed, and skip building a `rich` progress bar.

### Changed

- `populater --genome-server` no longer subscribes permanently. Its own help
  always said the URLs would not persist; they now really do not. Anyone who
  relied on the accidental persistence should run `refgenie subscribe -s <url>`
  once.
- The dashboard and the server are one application. `refgenie/server/main.py`
  is now a single factory, `create_app(mode="server"|"local")`, and
  `refgenie dash` runs it in local mode. The separate dash app
  (`refgenie.server.dash.main:app`) is gone.
- HTML is served by a React single-page app built into
  `refgenie/server/webui/`, replacing both Jinja template trees. A wheel ships
  the built bundle; a source checkout builds it with
  `npm --prefix frontend run build`. Until it exists, `/` answers 503 and the
  rest of the API is unaffected.
- The JSON API is served at `/v4` only. Both apps previously mounted it a second
  time at the root (`GET /genomes`, `GET /summary`, …); the root namespace now
  belongs to the web UI's client-side routes.
- `refgenie dash` binds `127.0.0.1` instead of `0.0.0.0`. Use `refgenie serve`
  for anything that has to be reachable from another machine.
- `GET /v4/assets/{digest}` and `GET /v4/assets` now include the asset's
  `seek_keys`.
- `GET /aliases`, `GET /aliases/{name}`, and `GET /assets/{asset_digest}/files`
  are now defined once, by the v4 server router. They were previously declared
  by both the v4 router and the core router with incompatible response schemas,
  so the published OpenAPI document advertised two contradictory definitions of
  each. The dash app no longer exposes those three paths.
- `GET /aliases` now returns a paginated response
  (`{"items": [...], "pagination": {...}}`) and honors `offset` / `limit` and
  the `q` / `search_fields` search parameters, matching every other collection
  endpoint. The previous unpaginated `{"items": [...]}` shape made the client's
  page loop repeat forever on a server with more aliases than one page.
- `psycopg` moved from the `dash` extra to a base dependency. `PostgresConfig`
  is a base configuration type, so a base install configured for PostgreSQL
  previously failed with a bare `ModuleNotFoundError`.

### Fixed

- Local mode now reads aliases from both the on-disk refget store and the SQL
  `alias` table. `refgenie store sync` writes federated aliases to SQL, so a
  build node that reads only its store cannot name any federated genome:
  `refgenie id <name>` failed for them and `refgenie catalog-export` published
  them as bare digests. The export now also refuses to run if it would drop
  alias rows the catalog holds.
- CLI flags that were accepted but never acted on are now either wired up or
  removed. Newly functional: `seekr -s/-p`, `listr -g/-s/-p`, `populater -p`,
  `remove --aliases`, `pull --force` on registry paths, `alias set --force`,
  `recipe requirements --recipe-version`, and
  `build --requirements --recipe-version`.
- `refgenie serve` on an install with `dash` but not `server` extras now prints
  the "server extras are not installed" message instead of raising
  `ModuleNotFoundError: No module named 'apscheduler'`.
- Building an asset from a thread other than the main thread no longer raises
  `ValueError: signal only works in main thread of the main interpreter`. The
  builder's SIGINT registration is now guarded the way the puller's already
  was; the CLI is unaffected.
- `build_asset(stage=True)` with no `genome_stage_folder` configured now raises
  before the build starts rather than after it finishes.
- The "replace the existing asset directory?" prompt in `pull` now goes through
  the caller's `confirm` callback instead of calling `rich.prompt.Confirm.ask`
  directly, so a non-interactive caller is refused instead of blocked on stdin.

### Removed

- The single-store `REFGENIE_REFGET_STORE_URL` environment variable and the
  `refget_store_url` argument to `Refgenie()` / `create_app()`. The server now
  federates over the `store` registry table instead; set `REFGENIE_STORES` to
  register stores on boot.
- `--genome` on `seekr`, `remove`, `rename`, `id` and `build`. Each of these
  commands takes registry paths, which already carry the genome.
- `--recipes` on `list` and `asset list`. Use `refgenie recipe list`.
- `--skip-asset-class` and `--skip-recipe` on `pull`. Asset classes and recipes
  arrive through data channels, not through `pull`.
- `--remote-class` on `seekr` and `populater`. Remote resolution always returns
  an https URL; `populater` versus `populate` is the only real distinction.
- `--text` on `build`. The `--requirements` table already degrades to plain text
  when stdout is not a terminal.
- `--requirements`, `--remote`, `--genome-server` and `--append-server` on
  `recipe show`; `--remote`, `--genome-server` and `--append-server` on
  `asset-class show`.
- The `recipe pull`, `recipe listr`, `recipe test`, `asset-class pull` and
  `asset-class listr` subcommands, which were never implemented.
- The `/v1/*` dash JSON API and the `/page/*` splash pages, with no redirects.
  The `/v1` prefix is reused for the local-only surface (`/v1/remote/*`, and the
  actions and jobs APIs that follow).
- The dash's private remote-catalog client and its TTL cache, including
  `GET /remote/assets`, `GET /remote/list`, `GET /remote/cache/status` and
  `POST /remote/cache/clear`. Remote browsing is `GET /v1/remote/genomes` and
  `GET /v1/remote/genomes/{digest}/assets`, which go through the same client the
  CLI uses. The hardcoded `http://refgenomes.databio.org/` fallback subscription
  is gone: no subscriptions now means an empty list.
- The `/static` mount. Brand assets ship with the web UI bundle.
- `DatabaseDialect` and `PostgresConfig.dialect`. The enum named `psycopg2`
  while the connection URL builds `postgresql+psycopg` (psycopg 3), and the
  field was never read.

## [1.0.0a1] - 2026-07-25

First public alpha of refgenie 1.0. This is a ground-up rewrite that supersedes
the legacy `refgenie` 0.x, `refgenconf`, and `refgenieserver` packages. It is
not compatible with 0.x configuration files or genome folders; see
[docs/migration-from-legacy.md](./docs/migration-from-legacy.md).

Highlights of the 1.0 design:

- A relational metadata store (SQLite or PostgreSQL) replaces the 0.x
  `genome_config.yaml`, with alembic migrations for schema changes.
- Content-addressed assets: an asset's digest covers the files a recipe declares
  it covers, and names are tracked separately in an `assetname` table.
- Sequence-collection (seqcol / refget) identifiers for genomes, with a refget
  store backend and GA4GH DRS endpoints on the server.
- Data channels: asset classes and recipes are synced from external indexes
  rather than shipped in the package.
- Separate serving modes — archive (tarball) and file-level — with staging and
  cloud push.
- An MCP server (`refgenie-mcp`, and `/mcp` on the HTTP server).

[Unreleased]: https://github.com/refgenie/refgenie/compare/v1.0.0a1...HEAD
[1.0.0a1]: https://github.com/refgenie/refgenie/releases/tag/v1.0.0a1
