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
[docs/development.md — "Filesystem side effects of the ORM"](./development.md#filesystem-side-effects-of-the-orm);
read that before writing code that deletes rows.

## Route ownership across app modes

There is one FastAPI factory, `refgenie/server/main.py::create_app(mode=...)`,
and it builds two apps:

| Prefix | Contents | Modes |
| --- | --- | --- |
| `/v4` | the shared JSON router; plus `version4` (archives, file downloads, summaries) in server mode | both |
| `/v1` | the local-only surface: `/v1/remote/*`, `/v1/jobs/*`, and `/v1/actions/*` (the state-changing command API; see `refgenie/server/local/`) | local |
| `/ga4gh/drs`, `/data_channel` | GA4GH and data-channel paths, at their root mount and under `/v4` | server |
| `/seqcol`, `/mcp` | mounted sub-applications | server |
| `/service-info` | GA4GH discovery, and the web UI's bootstrap document | both |
| `/` | the React SPA: its assets under `/_app/`, every other path falling back to `index.html` | both |

Two rules keep that table honest.

**One handler per path.** Every HTTP path is defined by exactly one handler in
exactly one router. FastAPI resolves collisions first-match-wins, so a path
defined twice is served by whichever router was included first, with the other
copy dead at runtime yet still advertised in the published OpenAPI document —
historically with two incompatible schemas. Paths both modes serve
(`GET /aliases`, `GET /aliases/{name}`, `GET /assets/{asset_digest}/files`, and
the entity listings) live on the shared router; server-only paths live on
`version4`.

**The root namespace belongs to the SPA.** Both routers used to be mounted a
second time at `""`, which put the JSON API in the same namespace as the UI's
client-side routes. Those mounts are gone: `GET /genomes` is the SPA's genome
page, and `GET /v4/genomes` is the API. The SPA catch-all is registered *last*
(it matches everything, so anything after it is dead) and answers an unmatched
path under an API prefix with a JSON 404 rather than an HTML 200.

`tests/test_web.py` asserts the mode isolation, the SPA-route/API
non-collision guard, and — with `tests/test_server.py` — that no route or
operationId is ever registered twice again.

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
