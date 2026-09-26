# Design notes

Short explanations of decisions that are not obvious from reading the code.

## Module naming: `models.py` vs `schemas.py`

Repeated module basenames are fine here and used on purpose — `manager.py`,
`const.py`, `helpers.py`, `queries.py` and a dozen others appear more than once,
scoped by their package. `models.py` is the exception, because the word is
ambiguous in a way the others are not: it can mean domain types, ORM tables, or
the JSON an HTTP client sends.

So one rule, enforced by `test_models_names_one_thing` in
`tests/test_layering.py`:

- **`refgenie/models.py`** is the only `models.py` in the package. It holds the
  domain types — `BuildParams`, `AssetRegistryPathComponents`, the validated
  string aliases — that the managers and the CLI both build on.
- **`schemas.py`** is the name for wire models: request bodies, response
  envelopes, SSE frames. There is one per web-facing package
  (`server/schemas.py`, `server/jobs/schemas.py`, `server/local/schemas.py`),
  and the repetition is the point — the name says what the file is.
- Internal records that are **not** on the wire stay with the code that owns
  them rather than moving to a models module. `Job` and `EventRing` live in
  `server/jobs/manager.py` because every field is guarded by the manager's lock
  and nothing else may touch them.

This drifted once already: a rename pass established `schemas.py`, and a later
commit adding the jobs and actions packages introduced two fresh `models.py`
files that nobody noticed. Hence the test.

## SQL table names

Following the suggestion in the SQLModel docs, table names are singular — `asset`,
`genome`, `alias`, and so on — to make the code more readable. See the
[SQLModel docs on SQL table names](https://sqlmodel.tiangolo.com/db-to-code/?h=singular##sql-table-names).

## Cascading deletes

In the database that holds refgenie asset metadata, deletes cascade to *child*
tables (`cascade="all, delete-orphan"` on the relevant relationships):

- `genome` -> `alias`, `asset_group`
- `asset_group` -> `asset`
- `asset` -> `staged_asset`, `seek_key`, `asset_name`
- `asset_class` -> its seek keys and links
- `recipe` -> its asset-class inputs
- `configuration` -> its owned rows

Deletes never cascade *upward*: removing every asset in a genome leaves the
`genome` row intact.

## Event listeners

On top of the cascading deletes, SQLAlchemy event listeners remove the
corresponding *files* on disk when rows are deleted. Deleting a genome therefore
deletes its assets and their files.

The listeners fire on flush but do not touch the disk there: they resolve the
paths — which is only possible while the row is still present — and queue them.
A session-level `after_commit` listener runs the queue; `after_rollback`
discards it. So `session.delete(obj); session.commit()` destroys data on disk
*after* the `COMMIT`, and a transaction that rolls back leaves every byte in
place. Every listener is enumerated in
[docs/refgenie/development.md — "Filesystem side effects of the ORM"](refgenie/development.md#filesystem-side-effects-of-the-orm);
read that before writing code that deletes rows.

## AssetManager is composed, not mixed

`AssetManager` inherits only `ResourceManager` and holds five smaller
managers, each taking its dependencies in its constructor. Do not split it with
mixins: a mixin reads host fields it never declares and calls back into host
methods, so it cannot stand alone, and the split follows the file rather than
the job. None of them holds a reference back
to `AssetManager`, and `AssetManager` holds no reference to the builder or the
puller:

- `rgc.asset.content` (`AssetContentManager`) -- the content write path:
  `add` (build output, pulled archive, or a folder the user registers),
  `adopt_name` and `add_incomplete`. It was the last large piece of
  `AssetManager`, and the only user of the asset-class manager. Its reads go
  straight to `queries.py`, inside its own session where one is open.
- `rgc.asset.group` (`AssetGroupManager`) -- asset groups and each group's
  default asset. Content, seek keys, the tree, the builder and the puller all
  depend on it; it depends on nothing asset-side.
- `rgc.asset.seek_key` (`SeekKeyManager`) -- seek-key lookups, and
  `resolve`, the one place a registry path's missing asset and seek key are
  filled in (local seek and remote seek's local shortcut both use it).
- `rgc.asset.tree` (`AliasTree`) -- the per-alias `alias/` and `builds/`
  trees on disk. `rgc.asset.content` renders the tree for every name it
  commits (`add`, `adopt_name`); only rename, set-alias, and the pull repair
  branch call `tree.render` / `render_genome` directly.
- `rgc.asset.links` (`AssetLinkManager`) -- parent/child links between assets.

Asset methods take the genome as a `genome_digest`; see "Genome identifiers
are typed" below.

**The builder and the puller do not belong to `AssetManager`.**
The build manager (`rgc.build`, `BuildManager`) and the puller (owned by
`rgc.transfer`, private `_puller`) are built lazily and depend on `rgc.asset`
one way. `AssetManager` does not construct them or hand them
`asset_manager=self`: that loop would let neither be understood alone. The build
manager is public and owns build (`run`), `preflight`, asset naming, the
Snakemake target templates and genome-init-plus-build
(`initialize_and_build`). It never names pull: the root hands it a
`pull_parent` callable for `run(pull_parents=True)`. Pulling is
`rgc.transfer.pull` / `pull_genomes` / `mirror` on the transfer manager.

**`Refgenie` has no mixins.** Database lifecycle is `rgc.database`, sequence
retrieval is `rgc.sequence`, pulling is `rgc.transfer`, building is
`rgc.build`. Each is a composed manager with explicit constructor
dependencies, not a mixin that silently assumes attributes of the root. The
Snakefile template and the web build form call `rgc.build` directly
(`target_template`, `default_asset_name`, `preflight`); there are no build
wrappers on the root. The build manager and the puller take `asset` and read
`asset.content`, `asset.group`, `asset.tree`, `asset.links` and the folders
from it. Data-channel sync is `rgc.sources.sync_channel`: it reads the channel
and registers into the asset-class and recipe managers it was handed.

**Build preparation is written once.** `BuildManager._prepare` fills recipe
defaults, resolves input assets, makes the `BuildCommandValues` and validates
the inputs. `BuildManager._build_command_values` makes the `BuildCommandValues`
for a build and for `preflight`. A build calls `_prepare` and stops at the
first problem; `preflight` calls the same steps one at a time, because it
reports each problem under its own field. `default_asset_name` is the one
naming path shared by `preflight` and the Snakefile; it uses the cached seek-key
probe. `_build` deliberately probes uncached, because what a build records has
to be probed at build time.

**`rgc.servers` is the pull side; "remote" means push targets.** `ServerManager`
owns the subscription list, one client per server, the read-only catalog
queries (`list_genomes`, `list_assets_for_genome`, `assets_table`) and remote
seek. It matches the `refgenie servers` command. It is also the single place
that asks servers what an alias or a collection is: `resolve_alias` (genome
source first, then the v4 client), `find_collection` and `genome_source`.
`set_genome_alias`, `GenomeManager.ensure_from_remote`, remote seek and
`refgenie id --remote` all call these rather than looping over servers
themselves. The word "remote" stays with
push targets (`refgenie remote`, the `Remote` table), which are `rgc.remote`
(see below). `rgc.sources` is data channels only. Creating a
genome from server metadata is a genome job, so it is
`GenomeManager.ensure_from_remote` / `init_from_remote`.

**The composed managers are public attributes, not wrapped.** Callers write
`rgc.asset.group.get_default(...)`, not `rgc.asset.get_default(...)` forwarding
to it. This is the same rule the `Refgenie` docstring states for the root: no
method whose only job is to delegate to one manager. A wrapper doubles the
surface that has to stay in sync and hides which object owns the behavior.

## Push remotes: `rgc.remote`

`RemoteManager` (`rgc.remote`) owns the `Remote` and `RemoteAssetLink` tables:
the places `refgenie push` uploads to, and which staged assets go to each.
Every reader and writer of the link table goes through its methods, never a
hand-written query, so the tests exercise the same path the product runs:
staging records push intent with `remote.link(...,
exist_ok=True)`, push reads `remote.unpushed` and calls `remote.mark_pushed`
after each verified upload, and the DRS and v4 routes get download URLs from
`remote.pushed_urls`. `publish_catalog.py` is the one exception: it copies raw
tables between databases to publish a catalog.

**One naming rule.** A remote is named by its id or its name, and nothing else:
an `int` or a string of digits is the id; any other string is the name. The
`name` column is required and unique, and `add` refuses a name that is empty,
all digits or already taken, so a name can never be mistaken for an id or match
two rows. Before this, `push --remote`
matched the name, `build --push-to` matched the prefix, and `remote remove`
matched the type, even though several remotes can share one.

## Genome identifiers are typed

A genome has two kinds of name: its digest (the sequence collection's level 0
digest, 32 characters) and any number of aliases (`hg38`). Both are strings, so
a bare `str` argument does not say which one it holds. Two small types in
`refgenie/models.py` do: `GenomeDigest` and `GenomeAlias`. They are `str`
subclasses with a Pydantic schema, and constructing one checks the value's
shape (a digest is 32 URL-safe base64 characters; an alias is non-empty and has
no `/`). The shape check only validates a value already known to be one kind.
Nothing uses it to decide which kind a string is.

The rule is that every function takes exactly one kind, stated by its type
hint:

- Library methods take `genome_digest: GenomeDigest` and use it as is.
- The one step from alias to digest is `AliasManager.resolve(GenomeAlias) ->
  GenomeDigest`, on every alias backend. It is alias-only: a string that is not
  an alias is a `MissingAliasError`, never a digest.
- Edges translate once. They know which kind they received and build that
  type: the CLI reads the genome in a registry path (and `-g`) as an alias and
  `--genome-digest` as a digest; web routes read `genome_name` and
  `genome_digest` fields; MCP tools take separate `alias` and `digest`
  parameters; `rgc.paths()` and the looper hook are keyed by alias, because
  pipeline configs name genomes that way. `build.run` takes a `GenomeAlias`,
  because the build folder is named after it.

There is no resolver that takes either kind and guesses (alias first, then "is
it a known digest?"), for two reasons. First, one string standing for two kinds
of value has no clean translation to a statically typed language, and a planned
Rust port needs the kinds apart. Second, nothing stops an alias from being
spelled like a digest, and a guess then silently picks the alias. With typed arguments that case
cannot be confused: a digest-shaped alias is only ever read as an alias, and a
digest only ever as a digest.

**`pull` is the one exception.** A pull must end with an alias, because the
alias tree it renders is keyed by one, but the genome may not exist locally
yet. An alias is used as given, or resolved on the server and registered; a
digest is registered under the aliases the server reports for it. So `pull`
takes `GenomeAlias | GenomeDigest` and branches on `isinstance`, the way a
Rust `enum GenomeRef { Alias, Digest }` would match. A plain `str` raises
`TypeError`, so the caller must say which it means. `transfer.pull_genomes`
follows it.

The types are notice for readers and callers, not enforcement: CI runs no type
checker, and nothing is wrapped in `@validate_call`. Values read back from the
database are plain `str`; they are wrapped where a digest enters a typed API
from outside (an edge, or a public return such as `AliasManager.resolve`), not
on every internal hop.

## Route ownership across app modes

There is one FastAPI factory, `refgenie/server/main.py::create_app(mode=...)`,
and it builds two apps:

| Prefix | Contents | Modes |
| --- | --- | --- |
| `/v4` | the `catalog` JSON router; plus `version4` (archives, file downloads, summaries) in server mode | both |
| `/v1` | the local-only surface: `/v1/remote/*`, `/v1/jobs/*`, and `/v1/actions/*` (the state-changing command API; see `refgenie/server/local/`) | local |
| `/ga4gh/drs`, `/data_channel` | GA4GH and data-channel paths, at their root mount and under `/v4` | server |
| `/seqcol`, `/mcp` | mounted sub-applications | server |
| `/service-info` | GA4GH discovery, and the web UI's bootstrap document | both |
| `/` | the React SPA: its assets under `/_app/`, every other path falling back to `index.html` | both |

Two rules keep that table honest.

**One handler per path.** Every HTTP path is defined by exactly one handler in
exactly one router. FastAPI resolves collisions first-match-wins, so a path
defined twice is served by whichever router was included first, with the other
copy dead at runtime yet still advertised in the published OpenAPI document,
possibly with an incompatible schema. Paths both modes serve
(`GET /aliases`, `GET /aliases/{name}`, `GET /assets/{asset_digest}/files`, and
the entity listings) live on `catalog`; server-only paths live on
`version4`.

**The root namespace belongs to the SPA.** No router is mounted at `""`: that
would put the JSON API in the same namespace as the UI's client-side routes.
`GET /genomes` is the SPA's genome page, and `GET /v4/genomes` is the API. The SPA catch-all is registered *last*
(it matches everything, so anything after it is dead) and answers an unmatched
path under an API prefix with a JSON 404 rather than an HTML 200.

`tests/server/test_app.py` asserts the mode isolation, the SPA-route/API
non-collision guard, and — with `tests/server/test_server.py` — that no route or
operationId is ever registered twice again.

## Server mode and the node's own store

`refgenie serve` always opens its database in server mode (`ServerMode` in
`refgenie/core/mode.py`): stores from the `store` registry, no sequence
ingestion. A database that was built on the same machine also has a
`<genome_folder>/.refget_store`, and that store is the only record of the
aliases and collections of the genomes built there. So when that folder exists,
server mode treats it as one more federated store, named `local`, ranked ahead
of every registered store, and reads aliases from it before the SQL table, the
same union local mode reads. A server with no store of its own reads the SQL
`alias` table alone.

The check needs the genome folder, which lives in the database, so the alias
manager is picked lazily on first use rather than in `Refgenie.__init__`. The
rejected alternatives: copying store aliases into SQL at build time or on serve
(two copies of one fact that drift), and choosing local mode whenever the
registry is empty (a build node that also registers remote stores would lose
its own names again).

## Background jobs in local mode

A pull takes minutes and a build takes hours, so neither can run inside a
request handler — and neither can run on FastAPI's shared threadpool, where one
build would occupy a slot for its whole duration and enough of them would
starve every synchronous JSON endpoint. `refgenie/server/jobs/` runs them on
two dedicated thread pools instead.

**Two pools, not one.** Builds are serialized in a single-slot executor because
pypiper's `PipelineManager` replaces `sys.stdout`/`sys.stderr` process-globally
for the duration of a pipeline; two overlapping builds corrupt that
save/restore chain. Pulls run two at a time: the `rich` live display that used
to forbid concurrency is not constructed when a `refgenie.progress` sink is
installed, and the job manager always installs one.

**How progress gets out.** `refgenie/progress.py` is a stdlib-only sink behind
a `ContextVar` — a no-op for the CLI and every other library consumer. Byte
counts come from `download_with_progress` through that sink; the managers'
narration ("Extracting asset tarball…") is picked up by a handler on the
`refgenie` logger and attributed to whichever job owns the calling thread; a
build's subprocess output is tailed off pypiper's own log file, because
pypiper's tee threads do not inherit the ContextVar.

**One stream, six event names.** `GET /v1/jobs/events` is a single SSE stream
multiplexed across every job, with a manager-global `seq` that `Last-Event-ID`
and `?since=` both resume from. There is no per-job stream: browsers cap
HTTP/1.1 at about six connections per origin and `EventSource` cannot set
headers. `GET /v1/jobs/{id}/events?format=json` is the polling fallback.

The six event names — `status`, `progress`, `log`, `done`, `truncated`,
`heartbeat` — are a closed set, not a convention. `EventSource` dispatches by
name and the client registers one listener per name it knows, so a seventh name
is not a compatible extension: the browser drops the frame before any JS sees
it. A coarse phase change is therefore a `progress` frame carrying a `phase`
and no byte counts, not a `stage` frame.

**A ref is not a record.** The 202 body of a submission is a `JobRef`, keyed
`job_id`; the state of a job is a `JobRecord`, keyed `id`. They are different
messages and the client has different types for them. `JobRecord` deliberately
does not inherit from `JobRef`, because that would rename one of the two.

**Phases are a shared vocabulary.** `JobProgress.phase` is a name from the list
in `frontend/src/components/jobs/phases.ts`, which turns it into a label and a
"step 4 of 8" counter — and, for the three phases that are multi-minute silent
hashing passes (`verify`, `digest`, `stage`), into explicit copy saying the
silence is normal. The backend emits a phase only where it can observe one
truthfully; an unrecognized or absent phase degrades to a plain label rather
than breaking. `tests/helpers.py::PHASE_VOCABULARY` mirrors the list so a typo
fails a test instead of silently losing the step counter.

**Jobs do not survive a restart.** They live in one process's memory, which is
also why local mode must run a single uvicorn worker. A user who restarts
`refgenie dash` mid-pull loses the job record; a partial `.tgz` may remain in
the asset group directory, and the next pull overwrites it.

## Integrations: the looper hook and Snakefile generation

`refgenie/integrations/` holds refgenie's side of outside workflow tools:
`looper.py` (the pre-submit hook) and `snakemake/` (Snakefile generation). These
are bridges, not front doors like `cli/`, `server/` and `mcp/`: each is reached
by a dotted path a user types into the other tool's config, or by one CLI
command. Nothing in `integrations/` imports looper or Snakemake, and
`tests/test_layering.py` holds it to that.

**The looper hook lives in refgenie.** Looper should not know about refgenie,
and a pipeline should not carry its own copy. Looper already calls any Python
function named in `pre_submit.python_functions`, so refgenie only needs to
provide one: `refgenie.integrations.looper.populate`. The looper hook is not a
refgenie plugin: the call goes from looper into refgenie. Calls the other way,
from refgenie out to other packages, go through the plugin system below.

**The hook is thin; the work is in `core/`.** `Refgenie.paths()`
(`core/paths.py`) is a general, lazy view of asset paths, and `Refgenie.populate`
(`core/populate.py`) resolves `refgenie://` strings. The hook wires the two into
looper's namespaces in a few lines.

**The path cache lives in memory only.** Looper calls the hook once per sample.
The hook builds one `Refgenie` and one path view per process and reuses them, so
each path is looked up at most once per run. There is no cache on disk: the
database already is the saved lookup, and a disk copy would go stale after a
pull or a remove. `paths()` is a method, not a cached property, so a
long-running process (the server) never serves paths from an old view.


## Plugins: refgenie calls out

`refgenie/plugins/` lets outside packages react to refgenie events. A plugin is
a function registered in an entry-point group named `refgenie.hooks.<hook>`.

**Five hooks.** `pre_pull` and `post_pull` fire once per asset pulled, bulk
pulls included. `pre_build` and `post_build` fire around a build. `post_update`
keeps its 0.x meaning, "local state changed": it fires once, after the
outermost public operation that changed local assets, aliases or genomes, and
only if something changed. It carries the list of changes. 0.x `pre_update`
fired before a YAML write that no longer exists, tags are gone, and running
plugins on `list` would slow every read, so those hooks were dropped.

**`hook(rg, event)`.** One frozen `HookEvent` rather than keyword arguments, so
refgenie can add fields without breaking a plugin. The return value is
ignored; no plugin can veto an operation.

**Managers record, the root dispatches.** Managers never load or call plugins.
A mutating manager method carries `@update_scope` and makes one
`record(Change(...))` call after its commit; the puller and the builder emit
their pre/post events. All of it goes into one `EventSink` that `Refgenie`
owns, and `PluginHost` (`plugins/host.py`) is the only code that calls a
plugin. Managers are marked rather than the root's methods because the CLI,
server and MCP call `rg.asset.remove`, `rg.asset.rename`,
`rg.asset.group.set_default`, `rg.genome.remove` and `rg.alias.remove`
directly. A manager built without a root gets `NULL_EVENTS` and does nothing.

**Why not SQLAlchemy events.** Local-mode aliases live in the RefgetStore, not
in SQL, so SQL events would miss alias changes. SQL events also fire on every
flush and commit, which would call `post_update` several times inside one pull.

**One `post_update` per operation.** `update_scope` nests: the outermost marked
call flushes, so a bulk pull of N assets gives N pre/post pull pairs and one
`post_update`, and `asset.content.add` fires after the alias tree is rendered. The
flush also runs when the operation raises, because committed changes are
real. Depth and pending changes are thread-local, since `refgenie dash` runs
jobs on worker threads that share one `Refgenie`. `PluginHost` has a
re-entrancy guard, so a plugin that pulls or builds does not fire hooks again.

**Layer 2.** `plugins/` imports only `config`, `logger` and `const`, so managers
(layer 5) and the root (layer 6) can both import it. Plugins get `rg` by duck
typing. `tests/test_layering.py` holds it there.

**Lazy discovery.** Nothing scans entry points at `import refgenie`, at
`Refgenie()`, or on a read. The first hook that fires makes one scan, cached
for the process, and loads only that hook's plugins.

**Log, never raise.** A plugin that fails to load or run is logged at WARNING
with its name and hook (traceback at DEBUG), and the operation and the other
plugins carry on. `Exception` is caught, not `BaseException`, so Ctrl-C works.

**The guard test.** `tests/plugins/test_update_scope_coverage.py` reads the
asset, group and genome managers statically: every public method that commits
must carry `@update_scope` or be named in its allowlist with a reason. That
keeps a new mutator from forgetting to fire `post_update`.

**Settings are one JSON column on `Configuration`.** refgenie1 has no
free-form config keys, so 0.x settings such as `nextflow_config:` had no
home, and one environment variable per plugin does not scale. Settings are
part of the configuration in force, and one column needs no joins, so they
are `configuration.plugin_settings`, a dict keyed by plugin name. It is on the
table class, not `ConfigurationPublic`, and the publish catalog strips it,
because settings may hold local paths. Changing a setting fires no hook.

**Server mode is off by default.** A deployed server serves many users, and a
plugin that rewrites a file there is almost never wanted, so
`Refgenie(server_mode=True)` runs no plugins unless `REFGENIE_SERVER_PLUGINS`
is set or the caller passes `plugins=True`. `refgenie dash` runs in local mode,
so plugins stay on there: a pull from the dash must update an nf-core config
exactly as `refgenie pull` does. `REFGENIE_DISABLE_PLUGINS` is the emergency
stop and wins over every opt-in.
