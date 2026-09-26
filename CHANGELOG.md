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

- `Refgenie.paths()`: a lazy, read-only view of local asset paths,
  `rg.paths()[genome][asset_group][seek_key] -> str`. Each path is looked up on
  first read and remembered for the life of the view; missing entries raise
  `KeyError`, so Jinja `is defined` guards work. `to_dict()` walks everything.
- Plugin system. Packages register functions in `refgenie.hooks.<hook>`
  entry-point groups (`pre_pull`, `post_pull`, `pre_build`, `post_build`,
  `post_update`). Each is called as `func(rg, event)`. `post_update` fires once
  after any operation that changes local assets, aliases or genomes. Failing
  plugins are logged, not raised. Plugins store settings in the database
  (`refgenie plugins set/unset`, `rg.plugins.settings(name)`; new
  `configuration.plugin_settings` column and migration). `refgenie plugins`
  lists plugins and settings. `REFGENIE_DISABLE_PLUGINS` turns plugins off.
  They are off by default in server mode (`REFGENIE_SERVER_PLUGINS=true` opts
  in).

### Changed

- Module moves, with no aliases kept: pulling lives in one package,
  `refgenie.managers.transfer` (`manager.py` for `TransferManager`, `puller.py`
  for `AssetPuller`, which moved from `refgenie.managers.asset.puller`, plus the
  pull-only prompt and SIGINT helpers). `refgenie.utils.build` is split into
  `refgenie.utils.digest` (content digests, build provenance, `checksum`) and
  `refgenie.utils.paths` (`get_build_dir`). `refgenie.catalog_transfer` is
  renamed `refgenie.publish_catalog`.

- Building is a public manager, `rgc.build` (`BuildManager`, the renamed
  `AssetBuilder`, now in `refgenie.managers.build`), and data-channel sync is
  on `rgc.sources`. `Refgenie` keeps no build or sync methods. Library renames,
  with no aliases kept:
  - `rgc.build_asset(...)` -> `rgc.build.run(...)`.
  - `rgc.preflight_build(...)` -> `rgc.build.preflight(...)`.
  - `rgc.initialize_and_build(...)` -> `rgc.build.initialize_and_build(...)`.
  - `rgc.get_asset_build_target_template(...)` -> `rgc.build.target_template(...)`.
  - `rgc.get_genome_init_target_template()` -> `rgc.build.genome_init_target_template()`.
  - `rgc.resolve_custom_seek_keys(recipe)` / `Refgenie.resolve_default_asset(...)`
    -> `rgc.build.resolve_custom_seek_keys(recipe)` / `BuildManager.resolve_default_asset(...)`.
    New: `rgc.build.default_asset_name(recipe)`.
  - `rgc.sync_data_channel(...)` -> `rgc.sources.sync_channel(...)`.
  - `refgenie.managers.asset.builder.AssetBuilder` -> `refgenie.managers.build.BuildManager`.
  A Snakefile generated before this change calls the old names; regenerate it.
- `Refgenie` has no mixins. Database lifecycle, sequence retrieval, and pulling
  are composed managers: `rgc.database` (`DatabaseManager`), `rgc.sequence`
  (`SequenceManager`), and `rgc.transfer` (`TransferManager`). CLI commands are
  unchanged. Library renames, with no aliases kept:
  - `Refgenie.get_database_config` -> `refgenie.managers.database.load_database_config`.
  - `Refgenie.get_default_database_engine` -> `refgenie.managers.database.create_database_engine`.
  - `rgc._create_db_and_tables()` -> `rgc.database.create_tables()`.
  - `rgc.init(...)` -> `rgc.database.init(...)` (returns `bool`); `rgc.init_backend` removed.
  - `rgc.check_for_db_migrations(log=)` -> `rgc.database.needs_migration(log=)`.
  - `rgc.migrate_db()` -> `rgc.database.migrate()`.
  - `rgc.check_table_exists` removed.
  - `rgc.purge(...)` -> `rgc.database.purge(...)`.
  - `rgc.getseq(genome_digest, locus)` -> `rgc.sequence.get(genome_digest, locus)`.
  - `refgenie.core.sequences._parse_locus` -> `refgenie.managers.sequence.parse_locus`.
  - `rgc.pull(...)` -> `rgc.transfer.pull(...)`.
  - `rgc.pull_all(genomes, ...)` and `rgc.pull_asset_for_genomes(asset_name, genomes, all_genomes, ...)`
    -> `rgc.transfer.pull_genomes(genomes, asset_group_name=..., all_genomes=...)`.
  - `rgc.init_genomes(...)` -> `rgc.transfer.init_genomes(...)`.
  - `rgc.mirror(...)` -> `rgc.transfer.mirror(...)`.
  - `refgenie/core/lifecycle.py`, `bulk.py`, and `sequences.py` are deleted.
- The asset content write path is its own manager, `rgc.asset.content`
  (`AssetContentManager`). Server-side alias and collection lookups all go
  through `rgc.servers`. Renames:
  - `rgc.asset.add_from_path` -> `rgc.asset.content.add`.
  - `rgc.asset.adopt_name` / `add_incomplete` -> `rgc.asset.content.adopt_name`
    / `add_incomplete`.
  - `rgc.check_remote_digest(digest, remote_servers=...)` ->
    `rgc.servers.find_collection(digest, server_urls=...)`.
  - New: `rgc.servers.resolve_alias(alias, server_urls=None)` asks the servers
    which genome an alias names, and `rgc.servers.genome_source(server_urls)`
    gives the first server usable as a genome source.
  - `rgc.asset.get(digest=...)` -> `rgc.asset.get_by_digest(digest)`, which
    returns None instead of raising. `rgc.asset.get` and `rgc.asset.exists`
    take a name only.
  - `rgc.populate_refgenie_registry_paths` and
    `rgc.populate_refgenie_registry_paths_in_string` are gone; use
    `rgc.populate` / `rgc.populater`.
  - `AssetGroupManager._remove` -> `remove_rows` (rows only, no disk cleanup;
    callers normally want `rgc.asset.remove`).
  - New: `rgc.asset.seek_key.resolve(...)` fills in the default asset and seek
    key and returns `(asset_name, seek_key)`.
  - `ServerManager`, `AssetBuilder` and `AssetPuller` take fewer constructor
    arguments: the builder and puller read `group`, `tree`, `links` and the
    folders from `asset`. `GenomeManager` requires its alias and server getters.
  - The eagerly built managers on `Refgenie` (`recipe`, `asset_class`,
    `configuration`, `remote`, `stage`, `sources`, `store`, `genome`) are plain
    attributes; the private `_recipe_manager`-style names are gone.
- Push remotes moved out of `ConfigurationManager` into `RemoteManager`
  (`rgc.remote`), which owns the `Remote` and `RemoteAssetLink` tables.
  Staging, `refgenie push` and the server's download redirects now all go
  through it. Renames:
  - `rgc.configuration.add_remote` / `upsert_remote` / `remove_remote` /
    `remote_table` / `remote_status` -> `rgc.remote.add` / `upsert` / `remove`
    / `table` / `status`. `add` takes `name=`, `type=`, `prefix=` (was
    `remote_type=`, `prefix=`, `description=`).
  - `rgc.configuration.link_asset_to_remote` / `unlink_asset_from_remote` /
    `mark_pushed` / `get_unpushed_links` -> `rgc.remote.link` / `unlink` /
    `mark_pushed` / `unpushed`. `link` returns `(link, created)` and takes
    `exist_ok=`; `unpushed` returns `(link, remote, staged_asset)` rows.
  - `refgenie.server.remote_assets.get_remote_url` -> `rgc.remote.pushed_urls`.
  - `rgc.configuration.remote_exists` is gone.
- A remote is named by its id or its name everywhere (`build --push-to`,
  `push --remote`, `remote status -r`, `remote remove`). Names are unique and
  may not be all digits.
- `refgenie remote add --description` is now `--name`.
- The `remote.description` column is now `remote.name`, required and unique.
- `refgenie remote remove` takes the remote's name or id as a positional
  argument instead of `--type`, which could not tell apart two remotes of one
  type.
- `refgenie push --remote` and `refgenie remote status -r` with a remote that
  does not exist now exit with "not found" instead of reporting nothing to do.
- The DRS remote access methods list only http/https remotes, each with its
  own URL. An s3 link used to be listed with the URL of whichever https remote
  also held the asset, and `/access/{type}:{id}:{mode}` ignored the id.

- Log messages and errors go to stderr, not stdout. Stdout carries only command
  output, so scripts that read `refgenie seek` output as a path no longer pick
  up an error message instead.
- `refgenie id --remote` looks up a digest given with `--genome-digest`
  (`refgenie id --genome-digest DIGEST --remote`). A positional name is always
  an alias; `id` no longer guesses from its shape that it is a digest.
  `--genome-digest` also works without `--remote`, alone for the genome or with
  an asset path (`refgenie id fasta --genome-digest DIGEST`).
- The store name `local` is reserved for the genome folder's own store, which
  server mode federates under that name. `refgenie store add local ...` and a
  `REFGENIE_STORES` entry named `local` are rejected.
- Genome identifiers are typed. Two `str` subclasses in `refgenie.models`
  (also importable from `refgenie`), `GenomeDigest` and `GenomeAlias`, say
  which kind of value a genome argument is. Nothing guesses any more whether a
  string is an alias or a digest, so an alias spelled like a digest can no
  longer shadow a genome.
  - Library methods take a digest, named `genome_digest=` (`genome_digests=`
    for lists). The `genome=` / `genomes=` argument that took either is gone,
    and so is `GenomeManager.resolve_digest`. The one way from an alias to a
    digest is `rgc.alias.resolve(GenomeAlias(...))`, which never falls back to
    treating its input as a digest. Old -> new:
    - `rgc.asset.get("fasta", "default", genome="hg38")` ->
      `rgc.asset.get("fasta", "default", genome_digest=rgc.alias.resolve(GenomeAlias("hg38")))`
    - `rgc.asset.group.get_default("fasta", genome=d)` ->
      `rgc.asset.group.get_default("fasta", genome_digest=d)`
    - `rgc.asset.table(genomes=[...])` -> `rgc.asset.table(genome_digests=[...])`
    - `rgc.getseq(genome=...)`, `rgc.asset.seek(genome=...)`,
      `rgc.servers.seek(genome=...)`, `rgc.add(genome=...)` -> `genome_digest=...`
    - `rgc.genome.resolve_digest(x)` -> `rgc.alias.resolve(GenomeAlias(x))`
      for an alias; a digest needs no resolving.
  - A registry path (`hg38/fasta`) names its genome by alias.
    `rgc.asset.seek_components` and `rgc.servers.seek_components` resolve it
    as one; `rgc.paths()` and the looper hook are keyed by alias.
  - `build_asset` and `preflight_build` take `genome_alias=` (was
    `genome_name=`), because the build folder is named after it.
  - `pull` is the one method that takes either kind, and says so in its type:
    `rgc.pull("fasta", genome=GenomeAlias("hg38"))` or
    `rgc.pull("fasta", genome=GenomeDigest(d))`. A plain `str` raises
    `TypeError`. `pull_all` and `pull_asset_for_genomes` take a list of either;
    `init_genomes` takes `aliases=` (was `genome_names=`).
  - The CLI treats the genome in a registry path, and `-g`, as an alias. A new
    `--genome-digest` flag names a genome by digest instead, which reaches a
    genome that has no alias: `refgenie seek fasta --genome-digest <digest>`.
    It is on `seek`, `seekr`, `getseq`, `list`, `listr`, `pull`, `add`,
    `remove`, `rename`, `stage add`, `stage remove`, `push`, `compare` and
    `genome remove`. Naming the genome both ways at once is an error.
    `genome set-metadata` keeps `--name` (alias) and `--digest`. `build` takes
    an alias only.
  - `DELETE /v1/actions/genomes/{genome_digest}` takes a digest only (the web
    UI already sent one). An alias or a malformed value there is a 404
    `genome_not_found`.
  - The MCP tools take explicit parameters instead of one string that could be
    either: `get_genome(alias=... | digest=...)`, likewise
    `get_genome_metadata`, `get_genome_sequences` and `list_assets`, and
    `compare_genomes(alias_a | digest_a, alias_b | digest_b)`. Exactly one of
    each pair is required.
- The Python API for assets and servers is reorganized; the old names are
  removed, not aliased. `rgc.asset` now holds asset records, the content write
  path and local seek, with four parts reached as attributes: `rgc.asset.group`,
  `rgc.asset.seek_key`, `rgc.asset.tree` and `rgc.asset.links`. Pulling from
  servers is `rgc.servers`; `rgc.sources` is data channels only. Old -> new:
  - `rgc.asset.get_group` / `group_exists` / `list_groups` / `get_default` /
    `set_default` -> `rgc.asset.group.get` / `exists` / `list_all` /
    `get_default` / `set_default`. `rgc.asset.remove_group` is gone from the
    public API (it left files behind); remove assets with `rgc.asset.remove`,
    which drops the group with its last asset.
  - `rgc.asset.get_seek_key` / `seek_key_exists` / `get_default_seek_key` /
    `list_seek_keys` -> `rgc.asset.seek_key.get` / `exists` / `get_default` /
    `list_all`.
  - `rgc.asset.render_alias_tree` -> `rgc.asset.tree.render_genome`;
    `rgc.asset.get_build_paths` / `find_build_dir` ->
    `rgc.asset.tree.build_paths` / `find_build_dir`. The `alias_whitelist=`
    argument is now `aliases=`.
  - `rgc.asset.get_parents` / `get_children` / `get_size` / `set_parents` ->
    `rgc.asset.links.parents` / `children` / `size` / `set_parents`.
  - `rgc.asset.seek_remote` -> `rgc.servers.seek`;
    `rgc.asset.list_remote` / `list_remote_assets_for_genome` /
    `list_remote_genomes` / `remote_table` -> `rgc.servers.list_assets` /
    `list_assets_for_genome` / `list_genomes` / `assets_table`;
    `rgc.asset.estimate_pull_size(...)` -> the function
    `refgenie.managers.sources.estimate_pull_size(...)`.
  - `rgc.asset.init_genome_from_remote` -> `rgc.genome.init_from_remote`.
  - `rgc.asset.build` / `pull` / `pull_multiple` are gone; use
    `rgc.build_asset`, `rgc.pull` and `rgc.pull_all` /
    `rgc.pull_asset_for_genomes` / `rgc.mirror`.
  - `rgc.configuration.subscribe` / `unsubscribe` / `get_server_subscriptions`
    and `rgc.sources.get_subscriptions` -> `rgc.servers.subscribe` /
    `unsubscribe` / `subscriptions`. `rgc.sources.get_server_client` ->
    `rgc.servers.client`; `rgc.server_clients` and `rgc.sources.server_clients`
    -> `rgc.servers.clients`. `rgc.sources.sync_server_clients` is removed.
  - For `refgenie://` paths already parsed into components:
    `rgc.asset.seek_components` (local) and `rgc.servers.seek_components`
    (remote).

  A bad client in `Refgenie(server_clients_mapping=...)` is now rejected when
  `rgc.servers` is first used, not in the constructor.
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
- The looper pre-submit hook is now `refgenie.integrations.looper.populate`
  (was `refgenie.looper_refgenie_populate_local`). Update the
  `pre_submit.python_functions` entry in your pipeline interface.
- The looper hook looks up asset paths lazily, once per run, instead of walking
  every genome, group and seek key for every sample.
- The looper hook no longer swallows errors while resolving `refgenie://` paths
  in the pipeline block. An unknown genome or a real failure now stops
  `looper run` with the actual error; a missing asset on a known genome still
  only warns and leaves the path as written (this now also holds for
  `refgenie populate`, which used to raise for a missing group).
- In the looper `refgenie` namespace, an asset group with no default asset is
  now undefined, matching `refgenie seek`. It used to fall back to any asset in
  the group.
- `refgenie.snakefile` moved to `refgenie.integrations.snakemake`.
- `psycopg` moved from the `dash` extra to a base dependency. `PostgresConfig`
  is a base configuration type, so a base install configured for PostgreSQL
  previously failed with a bare `ModuleNotFoundError`.

### Fixed

- `refgenie build --push-to <name>` records push intent again. It matched the
  remote's prefix, not its name, so a name logged "Remote not found" and the
  asset was silently never queued for push.

- `refgenie pull --init` exits non-zero when it registers no genome, instead of
  reporting success.
- `AssetManager.add_incomplete` and asset-group removal no longer open a second
  database session inside an open one. `add_incomplete` into an existing group
  no longer fails while logging the new asset.
- `refgenie serve` over a database built on the same machine lists and resolves
  its genomes' aliases again. Serving always runs in server mode, which read
  aliases only from the SQL `alias` table, while local builds keep them in the
  genome folder's `.refget_store`; clients could not pull a served genome by
  name. Server mode now federates over that store too, ahead of every
  registered store, and reads its aliases first.
- `refgenie pull --all` and `refgenie pull --asset` exit non-zero (2) when
  nothing could be resolved or pulled, instead of reporting success.
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
- `Refgenie.add` removed; use `rgc.asset.content.add`, which now renders the
  alias tree itself.

## [1.0.0a1] - 2026-07-25

First public alpha of refgenie 1.0. This is a ground-up rewrite that supersedes
the legacy `refgenie` 0.x, `refgenconf`, and `refgenieserver` packages. It is
not compatible with 0.x configuration files or genome folders; see
[docs/refgenie/migration-from-legacy.md](./docs/refgenie/migration-from-legacy.md).

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
