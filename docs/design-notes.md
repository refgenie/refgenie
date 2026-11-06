# Design notes

Short explanations of decisions that are not obvious from reading the code.

## Module naming: `models.py` vs `schemas.py`

Repeated module basenames are fine here (`manager.py`, `const.py`,
`helpers.py`, `schemas.py` and others appear more than once, scoped by their
package). `models.py` is the exception, because the word is ambiguous: it can
mean domain types, ORM tables, or the JSON an HTTP client sends.

So one rule, enforced by `test_models_names_one_thing` in
`tests/test_layering.py`:

- **`refgenie/models.py`** is the only `models.py` in the package. It holds the
  domain types (`BuildParams`, `AssetRegistryPathComponents`, the validated
  string aliases) that the managers and the CLI both build on.
- **`schemas.py`** is the name for wire models: request bodies, response
  envelopes, SSE frames, MCP tool results. There is one per web-facing package
  (`server/schemas.py`, `server/jobs/schemas.py`, `server/local/schemas.py`,
  `mcp/schemas.py`), so the name says what the file is.
- Internal records that are **not** on the wire stay with the code that owns
  them rather than moving to a models module. `Job` and `EventRing` live in
  `server/jobs/manager.py` because every field is guarded by the manager's lock
  and nothing else may touch them.

## SQL lives only in refgenie/db

Every SQL statement is written in `refgenie/db/`, and nowhere else. Services
(`server`, `cli`, `mcp`, `integrations`) call managers; managers call the db
layer. `refgenie/db/` holds the engine (`engine.py`), the session lifecycle
(`session.py`), schema creation and inspection (`schema.py`), the publish
catalog's table copying (`publish.py`), and `queries/`, one module per entity.

**The db-function contract.** A function in `refgenie/db/queries/` takes the
session first and never opens, commits or rolls one back: the manager that
called it owns the transaction. It returns rows, counts, or `(rows, total)`
for a page, never a statement. Every list has a fixed `ORDER BY` ending in a
unique column, so offset pages neither overlap nor skip. Each function
eager-loads everything its callers read, chosen by boolean keywords such as
`with_links=True` rather than loader objects, and rows stay usable after the
session closes (`expire_on_commit=False`). A database exception a manager must
react to is translated into a `refgenie.exceptions` type (`one_or_raise`,
`commit_or_conflict`).

**Managers** may hold a session to run a transaction (`add`, `delete`,
`flush`, `refresh`, `commit`, `rollback`), but never build or run a query.
They return table rows or small domain types from `refgenie/models.py`
(`Page`, `GenomeListing`).

**Services build their own wire shapes.** The HTTP API and the MCP tools call
the same manager reads (`rgc.genome.page`, `rgc.asset.get_by_digest`, ...)
and present them differently: `server/schemas.py` and the GA4GH page-based
envelope for HTTP, `mcp/schemas.py` and token paging with slim rows for MCP.
One query path, two presentations, so the two front doors cannot drift apart.

**The check.** `test_sql_is_written_only_in_the_db_layer` in
`tests/test_layering.py` fails on any import from `sqlalchemy` or `sqlmodel`
outside `refgenie/db/` that names anything but `Engine`, `Session`,
`SQLModel` or `Field`.

## SQL table names

Following the suggestion in the SQLModel docs, table names are singular (`asset`,
`genome`, `alias`, and so on) to make the code more readable. See the
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

The listeners fire on flush but do not touch the disk there. They resolve the
paths (possible only while the row is still present) and queue them.
A session-level `after_commit` listener runs the queue; `after_rollback`
discards it. So `session.delete(obj); session.commit()` destroys data on disk
*after* the `COMMIT`, and a transaction that rolls back leaves every byte in
place. Every listener is enumerated in
["Filesystem side effects of the ORM"](refgenie/development.md#filesystem-side-effects-of-the-orm);
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
  and `adopt_name`. It is the only user of the asset-class manager. Its reads go
  straight to `refgenie.db.queries`, inside its own session where one is open.
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
probe. `_build` probes uncached, because what a build records has to be probed
at build time.

**A recipe's image runs under Docker or Apptainer, through one module.**
`managers/containers.py` is the only code that knows how each engine runs an
image: `detect_engine` (Docker if its daemon answers, else Apptainer, or
`REFGENIE_CONTAINER_ENGINE`), `one_shot` for probes, `wrap` for build commands.
The recipe's `docker_image` is the single source of truth for where its tools
live, so a laptop (Docker) and a cluster with no Docker daemon (Apptainer) run
the same image and name the asset the same way. Each build command is its own
`docker run --rm` / `apptainer exec`, run by pypiper with `shell=True`: the
wrapped line quotes the recipe command, which may hold pipes that pypiper
would otherwise split.

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
two rows. Every command that takes a remote (`push --remote`, `build --push-to`,
`remote remove`) uses this rule.

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
of value has no clean translation to a statically typed language. Second, nothing stops an alias from being
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
| `/ga4gh/drs/v1` | GA4GH DRS, at the path the DRS spec requires (not under `/v4`) | server |
| `/data_channel` | data-channel paths, at their root mount and under `/v4` | server |
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

`tests/server/test_app.py` asserts the mode isolation and the SPA-route/API
non-collision guard, and with `tests/server/test_server.py` checks that no route
or operationId is registered twice.

## New configurations subscribe to the default server

`DatabaseManager.init` writes `DEFAULT_SERVER_URL` into the `servers` list of
the Configuration row it creates, so a first pull from `refgenie dash` (or over
the bridge, which may only pull from subscribed servers) works without setup.
Init records the URL and never contacts it, so it works offline. It
happens only when the row is created; re-running init keeps whatever the user
chose. Opt out with `refgenie init --no-default-server`, the
`subscribe_default_server=False` argument, or
`REFGENIE_SUBSCRIBE_DEFAULT_SERVER=false`. `tests/conftest.py` sets that
variable before importing refgenie, so no test (or CLI subprocess) can reach
the real server by accident; tests of the default pass the argument
explicitly.

Subscribing stays a same-origin action. The bridge can read the list
(`GET /v1/remote/servers`) but every "subscribe" affordance on the public site
is a link to the local dash's `/servers?subscribe=<url>`, which fills the form
and waits for a click.

## A pull brings its asset class and its parents

A fresh install has no asset classes, and a derived asset (`bwa_index`) links
to a parent (`fasta`) the client may not have. The puller fetches both from
the server it is pulling from, before downloading anything:

- **The class** is looked up by the server's `asset_class_name` (older servers:
  the group name) and, if no version of it is registered, fetched from
  `GET /v4/asset_classes/{id}` using the group's `asset_class_id`. That
  response carries the class's seek keys so the class can also be used for a
  local build. An older server sends no seek keys; the class is registered
  anyway with a warning, since the pulled asset carries its own seek keys
  from the server. A class registered this way is not rolled back if the pull
  fails: it is a definition, not data.
- **Missing parents** (by digest) are pulled first with `force=False` from the
  same server, each as a full pull of its own, so recursion and hooks come for
  free. `required_by` carries the chain of digests to stop a cycle. A parent
  that fails raises `ParentPullError` (a `PullFailedError`, so jobs report
  `pull_failed`) naming the parent and the reason. A parent that arrives is a
  complete asset and stays even if the child then fails; the child's
  transaction releases the parent's genome (`PullTransaction.release_genome`)
  so its rollback cannot remove the genome, and the parent with it.

In a bulk pull, an entry that is already present (often a parent an earlier
entry pulled) is logged and skipped rather than counted as a failure.

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
manager is picked lazily on first use rather than in `Refgenie.__init__`.
Copying store aliases into SQL would keep two copies of one fact that drift.
Choosing local mode whenever the registry is empty would lose a build node's
own names as soon as it registered a remote store.

A server that imports a publish catalog instead of building (api.refgenie.org)
has no store of its own, but the catalog still carries every genome built on
the build node, aliases included. Their collections live in the build node's
store, which the nightly mirrors to S3 as `refget-store/refgenie-main/`, so
that mirror must be in the server's `REFGENIE_STORES`, or the server lists
those genomes and answers 404 for their collections (and clients cannot pull
them). It is ranked after the curated stores, not ahead like `local`, so the curated
stores keep the default `/seqcol` mount. One mirror of the build node's store
is used instead of each built genome's registry store, because loading a whole
registry store for a few genomes costs a lot of memory, and each new store
would have to be added by hand.

## Background jobs in local mode

A pull takes minutes and a build takes hours, so neither can run inside a
request handler, nor on FastAPI's shared threadpool, where one
build would occupy a slot for its whole duration and enough of them would
starve every synchronous JSON endpoint. `refgenie/server/jobs/` runs them on
two dedicated thread pools instead.

**Two pools, not one.** Builds are serialized in a single-slot executor because
pypiper's `PipelineManager` replaces `sys.stdout`/`sys.stderr` process-globally
for the duration of a pipeline; two overlapping builds corrupt that
save/restore chain. Pulls run two at a time. That is safe because the `rich`
live display is not built when a `refgenie.progress` sink is installed, and the
job manager always installs one.

**How progress gets out.** `refgenie/progress.py` is a stdlib-only sink behind
a `ContextVar`, a no-op for the CLI and every other library consumer. Byte
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

The six event names (`status`, `progress`, `log`, `done`, `truncated`,
`heartbeat`) are a closed set. `EventSource` dispatches by
name and the client registers one listener per name it knows, so a seventh name
is not a compatible extension: the browser drops the frame before any JS sees
it. A coarse phase change is therefore a `progress` frame carrying a `phase`
and no byte counts, not a `stage` frame.

**A ref is not a record.** The 202 body of a submission is a `JobRef`, keyed
`job_id`; the state of a job is a `JobRecord`, keyed `id`. They are different
messages and the client has different types for them. `JobRecord` does not
inherit from `JobRef`, because that would rename one of the two.

**The request body is the job's parameters.** `POST /v1/actions/pull`, `/build`
and `/genomes` take `PullJobParams`, `BuildJobParams` and `GenomeInitJobParams`
as their bodies and hand them to the manager unchanged; the record's `params`
is the same object. Jobs are submitted only through `/v1/actions`.

**Phases are a shared vocabulary.** `JobProgress.phase` is a name from the list
in `frontend/src/components/jobs/phases.ts`, which turns it into a label and a
"step 4 of 8" counter. For the three long silent hashing phases (`verify`,
`digest`, `stage`) it also tells the user the silence is normal. The backend emits a phase only where it can observe one
truthfully; an unrecognized or absent phase degrades to a plain label rather
than breaking. `tests/support/actions.py::PHASE_VOCABULARY` mirrors the list so a typo
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
means "local state changed": it fires once, after the outermost public
operation that changed local assets, aliases or genomes, and only if something
changed. It carries the list of changes. There are no hooks on reads such as
`list`, because they would slow every read.

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

**Settings are one JSON column on `Configuration`.** Refgenie has no
free-form config keys, and one environment variable per plugin does not scale.
Settings are part of the configuration in force, and one column needs no
joins, so they are `configuration.plugin_settings`, a dict keyed by plugin
name. It is on the
table class, not `ConfigurationPublic`, and the publish catalog strips it,
because settings may hold local paths. Changing a setting fires no hook.

**Server mode is off by default.** A deployed server serves many users, and a
plugin that rewrites a file there is almost never wanted, so
`Refgenie(server_mode=True)` runs no plugins unless `REFGENIE_SERVER_PLUGINS`
is set or the caller passes `plugins=True`. `refgenie dash` runs in local mode,
so plugins stay on there: a pull from the dash must update an nf-core config
exactly as `refgenie pull` does. `REFGENIE_DISABLE_PLUGINS` is the emergency
stop and wins over every opt-in.

## Species names come from the stores

Species names are a flat table filled from the registry's `species.json`, not from NCBI. That needs no network and no lineage, and the file is the registry's own statement of what it calls each species.

Names are rows (`speciesname`), not a JSON list on the species. Name search is then one plain query that behaves the same on SQLite and Postgres, with no JSON-text matching or `\u` escaping.

There is one write rule everywhere: per-taxon replace. For each taxon written, its row is upserted and its names rewritten; other taxa are never touched. `store sync` merges every enabled store's file first, so the result does not depend on the order stores were synced in. A dropped name goes, and rows written by a pull from a server survive. A taxon that vanishes from every file keeps its rows; that is accepted, since a stale row only adds names to a taxon id.

The publish import rewrites species per taxon, the one exception to upsert-only. It is safe because no files hang off these tables, and Core deletes fire none of the ORM listeners that remove files. Upserting alone would keep a removed synonym matching on the server forever.

## Species names: search and display

Search: the genome `species` field matches the genome's own scientific name, plus every name of its taxon in the species list (an `IN` subquery on `speciesname`). A genome with no species row is still found by its scientific name.

Display: lists show one common name, the species list's. Other names are searchable everywhere but shown only in detail views (HTTP genome detail, MCP `get_genome`, the web UI genome page). HTTP list rows still carry the species object, because the web UI's species page groups and searches genomes in the browser.

The common name is looked up in `GenomeListing.common_name` at read time and never stored on the genome. A pull copies the server's species row, every time, so a pulled genome shows the same common name locally.

## Which FHR fields refgenie reads

`GenomeManager.apply_fhr_columns` copies these FHR fields onto the `genome` row, and nothing else:

| FHR field | genome column |
|---|---|
| `documentation` | `description` |
| `taxon.name` | `species_name` |
| `taxon.uri` | `taxon_id` |
| `genomeSynonym` | `synonyms` |
| `assemblySource` | `assembly_source` |
| `accessionID.name` | `assembly_accession` |
| `assemblyLevel` | `assembly_level` |

FHR `genome` is the genome's own name (`hg38`), not its species. It is not read, because the name already arrives as an alias (store collection aliases, `genome init --name`), and aliases are what every surface shows as a genome's name. Storing it again would only duplicate it. `version` (the assembly's version) is not read either. Both stay visible in the raw FHR record.

Common names come from the species list, not FHR. The spec has no common-name field and forbids extra ones, and the species list already carries every name of each species, keyed by taxon id. So there is no `genome.common_name` column; responses look the name up by `taxon_id`.

## Asset content digest

`get_dir_digest` hashes each file's relative path and bytes, because aligner indexes find their parts by filename. Symlinks (colocation links to parents), empty directories and file modes are left out. A recipe with no `inherent` list includes every file: including an extra file only splits a digest, which is visible, while leaving out a needed one makes different assets collide.

The build digest (`build_digest`) identifies the build rather than its output. Level 1 stores a digest of every attribute, not just the inherent ones, so the inherent set can change later by re-filtering stored rows.

## Progress reporting

`refgenie.progress` uses a `ContextVar` instead of a `progress=` parameter, so the sink reaches the download code without changing four public signatures. Each worker thread gets its own sink. Threads that library code starts do not inherit it, so build output is read by tailing pypiper's log file instead.

## Local-mode web security

`refgenie dash` listens on loopback, so any page the user visits can try to call it. `server/local/security.py` layers three guards: every action needs the `X-Refgenie-Action` header (which forces a CORS preflight against the origin allowlist); a `Host` guard accepts only loopback names (421 otherwise) to stop DNS rebinding; and CORS sends no credentials and never allows a cross-origin `DELETE`.

The public SPA may read a local instance (`bridge_mode=read`, the default) and, under `full`, call only `POST /v1/actions/pull` from a subscribed server. `/v1/files` is refused to bridge origins in every mode, so a public page can see which assets you have but never read them.

## A StagedAsset row is the publication record

Byte-serving and DRS routes check for a `StagedAsset` row (written only after staging finishes) and the asset's serving modes, not for files on disk, so a leftover tarball from an interrupted stage is never served.
