# Refgenie Test Suite

## Layout

```
tests/
├── conftest.py      # Shared fixtures and the tier-assignment hook
├── support/         # Shared plain functions, classes and constants
├── data/            # Test inputs (FASTA files, YAML configs)
├── scripts/         # services.sh, test-integration.sh, test-package.sh
├── e2e/             # e2e tier
├── integration/     # integration tier
└── <package>/       # one directory per refgenie/ subpackage (core, managers, db, ...)
```

Subdirectories mirror the `refgenie/` source tree, and each is a package, so a
basename may repeat across directories. Put a new test in the module that
already owns its subject rather than starting another file.

## Tiers

| Tier | Command | Needs |
|------|---------|-------|
| `unit` | `pytest` | nothing |
| `component` | `pytest -m component` | nothing |
| `e2e` | `pytest -m e2e` | `samtools` on PATH (some tests skip without it; see below) |
| `integration` | `./tests/scripts/test-integration.sh` | Docker, bulker |

No `samtools`? Use bulker to get it and run e2e locally with nothing skipped:
`bulker exec bulker_manifest.yaml -- python3 -m pytest -m e2e`. CI installs
`samtools` directly instead.

Ask pytest for counts and timings:

```bash
pytest --collect-only -q | tail -1        # how many tests in a selection
pytest --durations=10                     # where the time goes
```

Other useful selections:

```bash
pytest                                    # unit tier
pytest -m "unit or component"             # unit + component
pytest -m "unit or component or e2e"      # everything but integration
```

CI runs `unit`, `component` and `e2e` as separate jobs (see
`.github/workflows/test.yaml`). `integration` runs only from the script; run
through plain `pytest`, its tests skip.

### Which tier should I write in?

A test's tier comes from where it lives
(`pytest_collection_modifyitems` in `conftest.py`): anything under `tests/e2e/`
is `e2e`, anything under `tests/integration/` is `integration`, everything else
defaults to `unit`. Override per module or per class:

```python
pytestmark = pytest.mark.component  # whole module
```
```python
@pytest.mark.component  # single class
class TestCLIWithAssets: ...
```

Choose by what the test touches:

- **`unit`**: in-memory SQLite, `TestClient`, mocks, `tmp_path` scratch files.
  No subprocesses, external binaries, or network. This is the default and
  where new tests belong unless you have a reason.
- **`component`**: builds real genome folders, asset files or `.tgz` archives
  on disk and exercises them through the manager stack. Still no external
  binaries or network.
- **`e2e`**: shells out to the installed CLI (`subprocess`) or to real
  bioinformatics tools.
- **`integration`**: needs a live PostgreSQL or HTTP service.

If a test builds an asset, it is `component` at minimum. Tests that build
against a test double (see `managers/test_genome.py::TestGenomeInitFailure`) stay `unit`.

## Two homes: fixtures vs helpers

There is one rule about where shared test code lives, and the
hygiene tests in `tests/test_layering.py` enforce it:

- **`tests/support/`** — every plain importable function, class and constant,
  one module per concern (`factories`, `server`, `processes`, `remotes`,
  `snapshots`, `actions`). No fixtures. Its top-level imports are limited to the stdlib, `pytest` and
  `unittest.mock`; `refgenie`, `sqlmodel` and `fastapi` are imported inside the
  functions that need them, so the mock-only tests still run without the
  optional extras and `tests/e2e/cli_compat/` stays library-free.
- **`tests/conftest.py`** — fixtures and pytest hooks only. It imports from
  `tests/support/` what its fixtures need. The same rule holds for the
  nested conftests: test modules import plain helpers from the `helpers.py`
  beside them (`tests/integration/helpers.py`, `tests/e2e/cli_compat/helpers.py`),
  never from a conftest.

Never redefine a fixture name the root conftest already owns; the
hygiene tests in `test_layering.py` fail if you do.

## Key Fixtures

From `conftest.py`:

- **`engine`** - Function-scoped in-memory SQLite engine
- **`fixtures_path`** - Path to `tests/data/` (session-scoped)
- **`refgenie_minimal`** - Initialized Refgenie with no assets built;
  inits into the test's `tmp_path`, so its `genome_folder` differs per test
- **`refgenie_with_fasta`** - As above, plus the fasta asset class and recipe
- **`refgenie_fs`** - As above, plus the rCRSd genome initialized; nothing built
- **`refgenie_built`** - `refgenie_fs` with one fasta asset built as `test`
- **`staged_refgenie`** - `refgenie_built` plus staged records and a `.tgz`
- **`refgenie_session`** - Session-scoped Refgenie with built FASTA (read-only tests)
- **`server_client_world`** - `(client_rg, server_rg, url)`: a real server app
  served in-process and registered on the client's source manager
- **`bare_client_world`** - `(client_rg, server_rg, url, parent, child)`: like
  `server_client_world`, but the client has no asset classes (a fresh install)
  and the server holds an unstaged fasta plus an unstaged child linked to it
- **`fake_subprocess`** / **`failing_subprocess`** - `subprocess.run` patched to
  succeed silently, or (factory) to raise
- Autouse: `_isolated_default_genome_folder`, `reset_source_cache`

The `refgenie_*` fixtures that build assets are what separate `component` from
`unit`; if your test needs one of them, mark the module `component`.

## Shared helpers

From `tests/support/` (import them, do not re-implement them):

| Module | Helper | What it gives you |
|---|---|---|
| `factories` | `make_engine`, `register_fasta`, `build_rcrsd`, `stage_rcrsd` | catalog primitives |
| `factories` | `make_built_refgenie`, `make_server_rgc` | whole Refgenie worlds on disk |
| `server` | `make_server_app`, `make_local_app`, `make_server_client`, `make_local_client`, `serve_refgenie` | app + TestClient plumbing |
| `server` | `requires_server`, `requires_dash`, `needs_server` | the optional-extras guards (`needs_server` is the mark form, for a class or test) |
| `server` | `route_keys`, `effective_api_routes`, `assert_no_duplicate_routes`, `assert_unique_operation_ids` | FastAPI route inspection (the route-registration guards) |
| `processes` | `run_cli`, `popen_cli`, `cli_argv`, `find_free_port`, `wait_for_server` | the one CLI subprocess launcher |
| `factories`, `server` | `server_with_asset`, `seed_staged_asset`, `seed_asset`, `create_asset_class`, `create_asset_on_disk` | direct DB seeding for server tests |
| `server` | `make_server_client_world` | a server + client pair for pull/list-remote |
| `remotes` | `mock_server_client`, `real_client_no_init`, `MockRemoteSource`, `mocked_puller` | the pull/seekr mocks |
| `snapshots` | `fasta_asset`, `only_asset`, `content_dir`, `build_dir`, `build_flag`, `boom` | the standard built world |
| `snapshots` | `db_snapshot`, `disk_snapshot`, `assets_rows`, `asset_name_rows` | state snapshots |
| `factories` | `add_asset_from_files`, `stage_copy`, `make_command_values` | asset construction |
| `actions`, `factories` | `TESTS_DATA_DIR`, `RCRSD_FASTA`, `DEMO_FASTA`, `GENOME`, `GROUP`, `ASSET`, `OMIT`, `PUBLIC_ORIGIN`, `EVIL_ORIGIN`, `CAPABILITY_KEY_SET` | constants |

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
  coverage goes in the component tier (`tests/server/test_server.py`).
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
A bare `r.database.init()` falls back to `config.genome_folder`, which in production is
`~/.refgenie/genomes` -- real user data. Two safety nets exist and neither is a
license to omit the argument:

- `tests/conftest.py` sets `REFGENIE_HOME_PATH` to a `mkdtemp` directory *before*
  importing refgenie, so the process-wide default never points at `$HOME`.
- The autouse `_isolated_default_genome_folder` fixture repoints
  `config.genome_folder` / `config.genome_stage_folder` at the current test's
  `tmp_path`, so a bare `init()` (including the one
  `rgc.database.needs_migration` may issue on its own) lands in per-test scratch.

`tests/managers/test_database.py::test_bare_init_never_writes_under_home` guards both.

### Testing Server Endpoints

Use `make_server_client(rgc)`, which wraps `create_app(refgenie_instance=rgc)`
in a `TestClient`. Every router takes its dependencies through `Depends`, and
`create_app` overrides `get_refgenie` with the given instance -- there are no
globals to patch:

```python
from tests.support.server import make_server_client, server_with_asset


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
stops testing the one that ships, so a route can be broken in the real app
without a single failing test.

### Testing with Built Assets

These belong in the `component` tier. Use the session-scoped `refgenie_session`
where the test is read-only, so the FASTA is built once for the whole run:

```python
pytestmark = pytest.mark.component


def test_with_fasta(refgenie_session):
    # refgenie_session has the rCRSd genome with its fasta asset pre-built
    digest = refgenie_session.alias.resolve("rCRSd")
    path = refgenie_session.asset.seek(digest, "fasta")  # returns a str
    assert Path(path).exists()
```

The fasta recipe uses `refgenie-build-fasta` (pure Python), so no samtools is
needed. Reach for a function-scoped fixture only when the test mutates state.

## Refgenie Manager Pattern

The `Refgenie` class delegates to managers. Don't look for methods like `refgenie.get_aliases()` - use the manager attributes:

| Property | Manager | Common Methods |
|----------|---------|----------------|
| `refgenie.alias` | AliasManager | `resolve()`, `get_for_genome()`, `add()`, `remove()` |
| `refgenie.genome` | GenomeManager | `add()`, `get()`, `initialize_genome()` |
| `refgenie.asset` | AssetManager | `content.add()`, `get()`, `seek()`, `remove()` |
| `refgenie.servers` | ServerManager | `subscriptions()`, `subscribe()`, `client()`, `seek()`, `locate()`, `alias_names()` |
| `refgenie.sources` | SourceManager | `add_channel()`, `list_channels()` |
| `refgenie.remote` | RemoteManager | `add()`, `get()`, `link()`, `unpushed()`, `mark_pushed()`, `pushed_urls()` |

Example: `refgenie.alias.get_for_genome(genome_digest)` not `refgenie.get_aliases()`
