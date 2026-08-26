# Developer Quickstart

## Prerequisites

- **Python 3.11+**
- **[uv](https://docs.astral.sh/uv/)** - Python package manager
- **[Task](https://taskfile.dev/)** - Task runner

## Clone and install

```bash
git clone https://github.com/refgenie/refgenie.git
cd refgenie1
uv sync
```

To install with web UI extras:

```bash
uv sync --extra dash
```

To install with server extras (includes dash):

```bash
uv sync --extra server
```

## Running from source

Run any refgenie command through uv:

```bash
uv run refgenie <command>
```

For example:

```bash
uv run refgenie init
uv run refgenie list
uv run refgenie build genome_name/fasta --files fasta=genome.fa.gz
```


## Purge an existing local refgenie database

```bash
uv run refgenie purge -f
```

## Initialize a genome from local file

```bash
uv run refgenie init
uv run refgenie list  # no assets found
uv run refgenie asset-class add --source ../recipes/asset_classes/fasta_asset_class.yaml
uv run refgenie recipe add --source ../recipes/recipes/fasta_asset_recipe.yaml
uv run refgenie genome init --fasta tests/data/rCRSd.fa --name rCRSd -d "Human mitochondrial reference"
uv run refgenie build rCRSd/fasta --files fasta=tests/data/rCRSd.fa
```

Now, `uv run refgenie list` shows:

```sh
┏━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Aliases ┃ Genome digest                    ┃ Asset group ┃ Asset   ┃
┡━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━┩
│ rCRSd   │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ fasta       │ default │
└─────────┴──────────────────────────────────┴─────────────┴─────────┘
```

## Initializing a genome from a remote RefgetStore

Let's purge and re-initialize refgenie to reset, to show initializing from a remote store:


```bash
uv run refgenie purge -f
uv run refgenie init
```

You can initialize genomes from a remote RefgetStore (static files on S3). First, browse what's available:

```bash
uv run refgenie genome browse
uv run refgenie genome browse --server-url https://refgenie.s3.us-east-1.amazonaws.com/refget-store/jungle/
```

Then initialize a genome by its seqcol digest:

```bash
uv run refgenie genome init \
  --store https://refgenie.s3.us-east-1.amazonaws.com/refget-store/jungle/ \
  --digest="2Ls1P5eUdKbvtOhjJx3s2R5r0_I-IB5Z" \
  --name my_genome \
  -d "Description of the genome"
```

Note: if a digest starts with `-`, use `--digest="<DIGEST>"` (with `=`) so argparse doesn't treat it as a flag.


Once initialized, build an asset like this:

```bash
uv run refgenie data-channel add recipes http https://refgenie.github.io/refgenie-registry/index.yaml
uv run refgenie data-channel sync recipes --exists-ok
uv run refgenie asset-class add --source ../recipes/asset_classes/fasta_asset_class.yaml
uv run refgenie recipe add --source ../recipes/recipes/fasta_asset_recipe.yaml
uv run refgenie build my_genome/fasta
```



## Running the web UI in dev mode

The web UI is a React + Vite + TypeScript app in `frontend/`. Vite builds it
into `refgenie/server/webui/`, and the refgenie app serves it from there. A
source checkout has to build it once:

```bash
npm --prefix frontend install
npm --prefix frontend run build
```

Then start the local app:

```bash
uv run refgenie dash        # or: task dash
```

This binds `127.0.0.1:8080` and opens your browser. Without a built bundle the
app still starts and the API still answers; `/` returns a 503 page telling you
to run the build.

For UI work, run the two halves separately so both hot-reload:

```bash
# terminal 1 -- the API, with auto-reload
uv run uvicorn "refgenie.server.main:create_local_app" --factory --reload --port 8080

# terminal 2 -- vite, proxying /v4, /v1 and /service-info to port 8080
npm --prefix frontend run dev   # or: task dash-dev
```

Vite serves the UI on port 5173 and proxies API calls to the backend, so you
edit TypeScript and Python at the same time without rebuilding either.

To catch base-path and hashed-asset problems the way CI would, do a
production-style preview instead: `task web-build && uv run refgenie dash`,
then open `http://localhost:8080/`.

To ask a running server which UI build it is actually serving -- the built
assets are never committed, so this is the only way to know without checking
CI logs -- fetch `GET /build-info.json` directly, or read the `web_ui` block
under `refgenie` in `GET /service-info`. Both report the commit the bundle was
built from, whether that checkout was clean, and the build timestamp.


## Running MCP server

```bash
uv sync --extra server
uv run refgenie-mcp
```

## Running tests

### Unit tests

Fast tests using in-memory SQLite (~7s):

```bash
task test
```

Or directly:

```bash
uv run pytest -vv -s -x tests
```

### Integration tests

End-to-end tests requiring Docker PostgreSQL (~30s):

```bash
./tests/scripts/test-integration.sh
```

Do **not** run `pytest tests/integration/` directly -- the script handles starting and stopping Docker services.

### Running a single test

```bash
task ktest TARGET=alias
```

This runs all tests matching the pattern "alias".

## Filesystem side effects of the ORM

Refgenie keeps state in three places that must agree: the SQLite catalog, the
RefgetStore, and plain files on disk. Only the catalog has real transactions, so
it is the participant that decides.

**The rule: commit the catalog first, then clean up the filesystem, and make the
cleanup safe to re-run.** No code path may destroy managed data before the
catalog agrees it is gone. This applies to the event listeners below and equally
to any manager doing filesystem work of its own.

`refgenie/db/events.py` registers the mapper-level listeners (via
`register_events()`). The delete listeners run **on flush**, but they do not
touch the disk: they resolve the paths the delete makes obsolete — which is only
possible while the row is still present — and queue them through
`refgenie/db/cleanup.py`. Session-level listeners registered alongside them drain
that queue on `after_commit` and discard it on `after_rollback` /
`after_soft_rollback`.

So an ordinary

```python
session.delete(asset)
session.commit()
```

still deletes files and directories on disk, but only once the `COMMIT`
succeeds. Read this table before writing code that inserts or deletes these rows.

| Model | Event | What it does on disk |
| --- | --- | --- |
| `Asset` | `before_delete` | Queues `genome_folder / asset.path` — the asset's content directory — for removal after commit. Skipped when `asset.path is None` (incomplete asset); logs a warning and proceeds if the directory is already gone. |
| `Alias` | `before_delete` | Queues the trees `alias_owned_paths` names — `alias/<name>/` and `builds/<name>/` — for removal after commit. Only fires for the SQL-backed `AliasManager` (server mode); `StoreAliasManager` removes the same paths itself, after its own catalog (the store) has dropped the alias. |
| `StagedAsset` | `before_delete` | Queues, for after commit: `mode="archive"` the staged tarball; `mode="file"` the staging directory of per-file symlinks. Then prunes the parent group directory if it is empty. Originals in `genome_folder` are untouched. |
| `SeekKey` | `before_insert` | Read-only on disk. Stats the target file / directory / prefix glob to populate `seek_key.size`. Raises `ValueError` for an unsupported `SeekKeyType`. |
| *(session)* | `after_commit` | Runs the queued removals. Failures are logged, never raised: the catalog is already correct, so a file left behind is drift to sweep, not a reason to fail a command that succeeded. Must not emit SQL. |
| *(session)* | `after_rollback`, `after_soft_rollback` | Discards the queue. |

Two computations used to run in `before_insert` listeners and no longer do,
because both held the write lock for the length of a filesystem walk that had
nothing to do with the transaction:

- The staged tarball's `tarball_digest` / `tarball_size` are computed in
  `StageManager.create`, right after the tarball is written.
- An asset's `size` is computed in `AssetManager.add_from_path`
  (`utils.build.directory_size`), before the insert.

Consequences worth knowing:

- **A rollback no longer destroys data.** Rows come back and the files were never
  removed. The inverse — a committed delete whose file removal failed — leaves
  content with no catalog row, which is inert: no read path reaches it, and the
  write paths reclaim it (see `AssetManager._place_content`,
  `AssetPuller.pull`).
- **The path listeners resolve `genome_folder` from the newest `configuration`
  row**, not from the live `Refgenie` object. A test or script that points
  `genome_folder` somewhere unexpected will have these listeners act there.
- **Registration is global and idempotent.** `register_events()` guards on a
  module-level `_events_registered` flag, so listeners attach once per process
  and stay attached for every engine and session in it.

What this does *not* buy you: safety against concurrent writers. A second
process can still interleave with any of this. Single-process interruption is
recoverable; concurrency is a separate problem with a separate fix.

## Code formatting

```bash
task reformat
```

This runs [ruff](https://docs.astral.sh/ruff/) to format the codebase.

## Building a wheel

```bash
uv build
```

The wheel is created in the `dist/` directory.

To validate the built package includes all required files:

```bash
./tests/scripts/test-package.sh
```

This builds the wheel, installs it in an isolated environment, and runs smoke tests.

## Task runner commands

All available task runner targets:

| Command | Description |
|---------|-------------|
| `task test` | Run unit tests |
| `task ktest TARGET=<pattern>` | Run tests matching a pattern |
| `task reformat` | Format code with ruff |
| `task init` | Initialize refgenie backend |
| `task purge` | Purge backend and remove all assets |
| `task reinit` | Purge then re-initialize |
| `task dash` | Run the local web UI (serves the built bundle) |
| `task dash-dev` | Run the vite dev server for the web UI |
| `task web-install` | Install frontend dependencies (`npm ci`) |
| `task web-build` | Build the web UI bundle into `refgenie/server/webui/` |
| `task web-dev` | Run the vite dev server for the web UI (alias of `dash-dev`) |
| `task web-check` | Typecheck and lint the frontend |
| `task generate-snakefile` | Generate snakemake pipeline |
| `task archive JOBS_COUNT=<n>` | Run snakemake archive pipeline |
| `task alembic-revision MESSAGE=<msg> DB_CONN_STR=<url>` | Create new DB migration |
| `task alembic-upgrade DB_CONN_STR=<url>` | Upgrade DB schema to latest |
| `task build` | Build package distribution |
| `task package` | Build the wheel/sdist and validate them in an isolated venv |
