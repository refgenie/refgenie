# Refgenie Test Suite

## Directory Structure

```
tests/
├── conftest.py              # Shared fixtures + the tier-assignment hook
├── helpers.py               # Shared plain functions, classes and constants
├── data/                    # Test fixtures (FASTA files, YAML configs)
├── scripts/                 # Test runner scripts
│   ├── services.sh          # Docker PostgreSQL and HTTP server management
│   └── test-integration.sh  # Integration test runner
├── e2e/                     # e2e tier
│   └── cli_compat/          # Python-vs-Rust CLI compatibility, all subprocess
├── integration/             # integration tier (Docker, external services)
│   ├── conftest.py          # Integration-specific fixtures + server subprocesses
│   ├── test_cli_integration.py
│   ├── test_library_integration.py
│   ├── test_pull_integration.py
│   └── test_server_integration.py
└── test_*.py                # unit tier, plus the modules marked `component`
```

### Module inventory

One subject per module. Put a new test in the module that owns its subject
rather than starting another file.

| Module | Subject it owns | Tier |
|---|---|---|
| `test_asset.py` | Asset CRUD/registry-path parsing, read-only queries, asset-class lifecycle, seek-key builder/CLI/persistence, `refgenie://` population and `add`/insert, non-path seek keys and seekr file mode | mixed |
| `test_asset_content.py` | `add_from_path` write ordering, staging, colocation symlinks, digest-addressed asset names | component |
| `test_asset_removal.py` | Asset and genome removal ordering | component |
| `test_build_dir_location.py` | Pre-build input validation and build-directory location/flag invariants | component |
| `test_cli.py` | argv parsing, model validation, help rendering, exit codes, dispatch, `refgenie id`, the `refgenie-build-fasta` console script | mixed |
| `test_data_channels.py` | Data-channel CRUD and index-file parsing | unit |
| `test_db.py` | ORM round-trips, Alembic chain, manager facades, database config and migration state, error hierarchy, fixtures-vs-helpers hygiene | mixed |
| `test_dir_digest.py` | `get_dir_digest`: the asset primary key | unit |
| `test_genome.py` | GenomeManager, genome init/build flow, `get_metadata`, `Refgenie.init` invariants, alias resolve/add/remove and alias trees, `getseq` against a real RefgetStore | mixed |
| `test_lint.py` | Package hygiene (ruff over the tree) | unit |
| `test_mcp.py` | MCP tools and the `mcp`-extras import guard | unit |
| `test_pull.py` | `Refgenie.pull` end to end, rollback, `list_remote` | mixed |
| `test_recipes.py` | Recipe CRUD and validation, multi-version asset classes and recipes, Snakefile generation | mixed |
| `test_remote_push.py` | RemoteAssetLink, remote status, the push workflow, catalog export/import | mixed |
| `test_server.py` | The server HTTP app: endpoints, service info, data-channel router, route guards | mixed |
| `test_server_client.py` | `RefgenieserverClient`, store-mode selection, download progress | unit |
| `test_serving_modes.py` | Serving-mode resolution and staging by mode, GA4GH DRS endpoints vs serving modes | mixed |
| `test_sources.py` | Genome sources: service-info bootstrap, source cache, rgstore detection, remote genomes, browse/sync | unit |
| `test_utils.py` | `refgenie.utils.*`: IO/YAML, checksums, encryption, tarballs, progress columns, confirmation prompts | mixed |
| `test_web.py` | The whole web UI surface: the `refgenie.progress` hook, JobManager, pull/build error mapping, the jobs HTTP API and its multiplexed SSE stream, the local actions API `/v1/actions/*` and its security stack, the localhost bridge, `create_app(mode)` and mode isolation, SPA serving, the OpenAPI-to-TypeScript drift guard, threaded builds, and a real end-to-end pull/build | mixed |

## Tiers

| Tier | Command | Needs |
|------|---------|-------|
| `unit` | `pytest` | nothing |
| `component` | `pytest -m component` | nothing |
| `e2e` | `pytest -m e2e` | `samtools` on PATH (some tests skip without it) |
| `integration` | `./tests/scripts/test-integration.sh` | Docker, bulker |

Roughly: the unit tier is the inner loop and takes seconds, `e2e` and
`integration` take about a minute each, and everything together takes a few
minutes. Counts and wall times are deliberately not recorded here — they rot.
Ask pytest instead:

```bash
pytest --collect-only -q | tail -1        # how many tests in a selection
pytest --durations=10                     # where the time actually goes
```

Other useful selections:

```bash
pytest                                    # unit tier -- the inner loop
pytest -m "unit or component"             # the old default scope
pytest -m "unit or component or e2e"      # everything but integration
pytest -m "unit or component or e2e or integration"   # collect all four tiers
```

CI runs `unit`, `component` and `e2e` as separate jobs (see
`.github/workflows/test.yaml`). `integration` runs only from the script.

Never run `pytest tests/integration/` directly — tests will skip without the
setup script.

### Which tier should I write in?

A test's tier comes from where it lives
(`pytest_collection_modifyitems` in `conftest.py`): anything under `tests/e2e/`
is `e2e`, anything under `tests/integration/` is `integration`, everything else
defaults to `unit`. Override per module or per class:

```python
pytestmark = pytest.mark.component      # whole module
```
```python
@pytest.mark.component                  # single class
class TestCLIWithAssets: ...
```

Choose by what the test actually touches:

- **`unit`** — in-memory SQLite, `TestClient`, mocks, `tmp_path` scratch files.
  No subprocesses, no external binaries, no network. Should be well under 0.1s.
  This is the default and where new tests belong unless you have a reason.
- **`component`** — builds real genome folders, asset files or `.tgz` archives
  on disk and exercises them through the manager stack. Typically 0.2-0.5s each.
  Still no external binaries and no network.
- **`e2e`** — shells out to the installed CLI (`subprocess`) or to real
  bioinformatics tools. ~2.5s each, dominated by interpreter startup.
- **`integration`** — needs a live PostgreSQL or HTTP service.

If a test builds an asset, it is `component` at minimum. Tests that build
against a test double (see `test_genome.py::TestGenomeInitFailure`) stay `unit`.

## A note on `addopts = "-s"`

`pyproject.toml` passes `-s` (capture disabled). Historically the suite errored
en masse without it; that is no longer true — the full 1072-test scope passes
with `--capture=fd` and zero errors, and `-s` costs nothing measurable
(169.5s vs 168.5s over the full scope). The mechanism behind the old breakage
was pypiper's `PipelineManager`, which swaps `sys.stdout`/`sys.stderr` for
`_LogTee` wrappers at construction and restores them at `stop_pipeline()`;
older pypiper did this at the file-descriptor level, which collided head-on
with pytest's `--capture=fd`. pypiper 0.15.1 only touches the Python-level
streams, so the collision is gone. `-s` is kept for now because the noisy log
output is occasionally useful, not because it is required.

Related, still true: `PipelineManager.__init__` calls
`atexit.register(self._exit_handler)` and never unregisters. Every
asset-building test leaves one behind (measured: 37 queued after 30 tests), so
a full `component` run fires a few hundred of them at interpreter shutdown.
Harmless (~0.5s) but it explains the odd activity after the last test.

## Two homes: fixtures vs helpers

There is one rule about where shared test code lives, and the
fixture-hygiene tests in `tests/test_db.py` enforce it:

- **`tests/helpers.py`** — every plain importable function, class and constant.
  No fixtures. Its top-level imports are limited to the stdlib, `pytest` and
  `unittest.mock`; `refgenie`, `sqlmodel` and `fastapi` are imported inside the
  functions that need them, so the mock-only tests still run without the
  optional extras and `tests/e2e/cli_compat/` stays library-free.
- **`tests/conftest.py`** — fixtures and pytest hooks only. It imports from
  `tests/helpers.py` what its fixtures need.

Never redefine a fixture name the root conftest already owns. A local
`refgenie_built` once silently changed the asset name for ten tests, and a
duplicate `fixtures_path` shadowed the root one for the whole integration
suite; the fixture-hygiene tests in `test_db.py` now fail on either.

## Key Fixtures

From `conftest.py`:

- **`engine`** - Function-scoped in-memory SQLite engine
- **`fixtures_path`** - Path to `tests/data/` (session-scoped)
- **`refgenie_minimal`** - Initialized Refgenie with no assets built (~10ms);
  inits into the test's `tmp_path`, so its `genome_folder` differs per test
- **`refgenie_with_fasta`** - As above, plus the fasta asset class and recipe (~21ms)
- **`refgenie_fs`** - As above, plus the rCRSd genome initialized; nothing built
- **`refgenie_built`** - `refgenie_fs` with one fasta asset built as `test`
- **`staged_refgenie`** - `refgenie_built` plus staged records and a `.tgz`
- **`refgenie_session`** - Session-scoped Refgenie with built FASTA (read-only tests)
- **`server_client_world`** - `(client_rg, server_rg, url)`: a real server app
  served in-process and registered on the client's source manager
- **`fake_subprocess`** / **`failing_subprocess`** - `subprocess.run` patched to
  succeed silently, or (factory) to raise
- Autouse: `_isolated_default_genome_folder`, `reset_source_cache`

The `refgenie_*` fixtures that build assets are what separate `component` from
`unit`; if your test needs one of them, mark the module `component`.

## Shared helpers

From `tests/helpers.py` (import them, do not re-implement them):

| Helper | What it gives you |
|---|---|
| `make_engine`, `register_fasta`, `build_rcrsd`, `stage_rcrsd` | catalog primitives |
| `make_built_refgenie`, `make_server_rgc` | whole Refgenie worlds on disk |
| `make_server_app`, `make_local_app`, `make_server_client`, `make_local_client`, `serve_refgenie` | app + TestClient plumbing |
| `requires_server`, `requires_dash` | the two optional-extras guards |
| `run_cli`, `popen_cli`, `cli_argv`, `find_free_port`, `wait_for_server` | the one CLI subprocess launcher |
| `server_with_asset`, `seed_staged_asset`, `seed_asset`, `create_asset_class`, `create_asset_on_disk` | direct DB seeding for server tests |
| `make_server_client_world` | a server + client pair for pull/list-remote |
| `mock_server_client`, `real_client_no_init`, `MockRemoteSource`, `mocked_puller` | the pull/seekr mocks |
| `fasta_asset`, `only_asset`, `content_dir`, `build_dir`, `build_flag`, `boom` | the standard built world |
| `db_snapshot`, `disk_snapshot`, `assets_rows`, `asset_name_rows` | state snapshots |
| `add_asset_from_files`, `stage_copy`, `make_command_values` | asset construction |
| `TESTS_DATA_DIR`, `RCRSD_FASTA`, `DEMO_FASTA`, `GENOME`, `GROUP`, `ASSET`, `OMIT` | constants |

### Integration Test Fixture Conventions (`tests/integration/`)

- Expensive fixtures (genome builds, server startups) should be session-scoped
  and shared across test classes. Avoid duplicating a genome build when an
  existing fixture already covers it.
- Server fixtures live in `tests/integration/conftest.py` and serve multiple
  test classes (`multi_remote_servers`, `file_mode_pull_server`,
  `refgenie_serve_subprocess`, `store_backed_server`). Add new genomes/assets to
  an existing fixture via `build_server_env` rather than starting another server.
- Servers are real `refgenie serve` subprocesses over real sockets. Hand-rolled
  route fakes and in-process test apps do not belong in this tier; route-level
  coverage goes in the component tier (`tests/test_server.py`).
- Test classes should be stateless consumers of shared fixtures. Use
  function-scoped `tmp_path` only for test-specific scratch data (e.g.,
  extracting tarballs).

## Testing Patterns

### Testing Refgenie Operations

Use `refgenie_minimal` fixture for lightweight tests:

```python
def test_alias_operations(refgenie_minimal):
    rgc = refgenie_minimal
    rgc.genome.add("digest123", "Test genome", ["alias1"])
    assert rgc.alias.resolve("alias1") == "digest123"
```

**Every `init()` in a test passes an explicit `genome_folder` under `tmp_path`.**
A bare `r.init()` falls back to `config.genome_folder`, which in production is
`~/.refgenie/genomes` -- real user data. Two safety nets exist and neither is a
license to omit the argument:

- `tests/conftest.py` sets `REFGENIE_HOME_PATH` to a `mkdtemp` directory *before*
  importing refgenie, so the process-wide default never points at `$HOME`.
- The autouse `_isolated_default_genome_folder` fixture repoints
  `config.genome_folder` / `config.genome_stage_folder` at the current test's
  `tmp_path`, so a bare `init()` (including the one
  `check_for_db_migrations` may issue on its own) lands in per-test scratch.

`tests/test_genome.py::test_bare_init_never_writes_under_home` guards both.

### Testing Server Endpoints

Use `make_server_client(rgc)`, which wraps `create_app(refgenie_instance=rgc)`
in a `TestClient`. Every router takes its dependencies through `Depends`, and
`create_app` overrides `get_refgenie` with the given instance -- there are no
globals to patch:

```python
from tests.helpers import make_server_client, server_with_asset

@pytest.fixture
def server_client(tmp_path):
    rgc, digests = server_with_asset(tmp_path, assets=[...])
    with make_server_client(rgc) as client:
        yield client, digests

def test_endpoint(server_client):
    client, _ = server_client
    assert client.get("/v4/genomes/...").status_code == 200
```

For the local app (`refgenie dash`) use `make_local_client(rgc)`, which builds
the same `create_app` factory with `mode="local"`. Never hand-build a
`FastAPI()` out of routers to get a client: a test that assembles its own app
stops testing the one that ships, which is how the dash's entire browse path
was a hard 500 without a single failing test.

### Testing with Built Assets

These belong in the `component` tier. Use the session-scoped `refgenie_session`
where the test is read-only, so the FASTA is built once for the whole run:

```python
pytestmark = pytest.mark.component

def test_with_fasta(refgenie_session):
    # refgenie_session has the rCRSd genome with its fasta asset pre-built
    path = refgenie_session.asset.seek("rCRSd", "fasta")
    assert path.exists()
```

The fasta recipe uses `refgenie-build-fasta` (pure Python), so no samtools is
needed. Reach for a function-scoped fixture only when the test mutates state.

## Refgenie Manager Pattern

The `Refgenie` class delegates to managers. Don't look for methods like `refgenie.get_aliases()` - use manager properties:

| Property | Manager | Common Methods |
|----------|---------|----------------|
| `refgenie.alias` | AliasManager | `resolve()`, `get_for_genome()`, `add()`, `remove()` |
| `refgenie.genome` | GenomeManager | `add()`, `get()`, `initialize_genome()` |
| `refgenie.asset` | AssetManager | `add()`, `get()`, `build()`, `pull()`, `seek()` |
| `refgenie.sources` | SourceManager | `get_subscriptions()`, `add_channel()` |

Example: `refgenie.alias.get_for_genome(genome_digest)` not `refgenie.get_aliases()`
