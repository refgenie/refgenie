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
│   ├── helpers.py           # Plain helpers the integration tests import
│   ├── test_cli_integration.py
│   ├── test_library_integration.py
│   ├── test_pull_integration.py
│   └── test_server_integration.py
├── core/                    # refgenie/core/: store router, federation, removal
├── managers/                # refgenie/managers/: asset, genome, sources, recipes, pull, serving, remote
├── db/                      # refgenie/db/: tables, migrations, catalog config
├── utils/                   # refgenie/utils/: helpers, dir digest, build provenance
├── cli/                     # refgenie/cli/: parsing, dispatch, remote push
├── mcp/                     # refgenie/mcp/: MCP tools
├── plugins/                 # refgenie/plugins/: hook events, dispatch, registry (fake_plugin.py is a fixture plugin)
├── integrations/            # refgenie/integrations/: looper and Snakemake bridges
├── server/                  # refgenie/server/: HTTP app, app factory, jobs, local actions, bridge
└── test_*.py                # package-wide hygiene (layering, lint, config layout) and root modules (progress, publish catalog)
```

Subdirectories mirror the `refgenie/` source tree. Each is a package (it has an
`__init__.py`), so a basename may repeat across directories. A module's
directory does not set its tier (see "Which tier should I write in?").

### Module inventory

One subject per module. Put a new test in the module that owns its subject
rather than starting another file.

| Module | Subject it owns | Tier |
|---|---|---|
| `test_config_layout.py` | `refgenie/config/__init__.py` only indexes; definitions live in named modules | unit |
| `test_layering.py` | The package's dependency direction: imports point down or sideways, never up | unit |
| `test_lint.py` | Package hygiene (ruff over the tree) | unit |
| `test_progress.py` | The `refgenie.progress` hook: no-op without a sink, thread-scoped sinks, abort, the download progress sink | unit |
| `core/test_asset_removal.py` | Asset and genome removal ordering | component |
| `core/test_store_federation.py` | Federated multi-store serving: store registry, alias policy, sync, serving a locally built database | component |
| `core/test_store_router.py` | The federated `RefgetStoreRouter` | unit |
| `managers/test_asset.py` | Asset CRUD/registry-path parsing, read-only queries, asset-class lifecycle, seek-key builder/CLI/persistence, `refgenie://` population and `add`/insert, non-path seek keys | mixed |
| `managers/test_asset_alias_tree.py` | `AliasTree` (`rgc.asset.tree`): rendering the per-alias view, default-asset rendering, alias whitelists | component |
| `managers/test_asset_content.py` | `AssetContentManager` (`rgc.asset.content`): `add` write ordering, staging, colocation symlinks, digest-addressed asset names | component |
| `managers/test_asset_group.py` | `AssetGroupManager` (`rgc.asset.group`): group lookups, the one-default-per-group invariant | mixed |
| `managers/test_asset_seek_key.py` | `SeekKeyManager` (`rgc.asset.seek_key`): seek-key lookups and the default seek key | unit |
| `managers/test_build_dir_location.py` | Pre-build input validation and build-directory location/flag invariants | component |
| `managers/test_data_channels.py` | Data-channel CRUD and index-file parsing | unit |
| `managers/test_genome.py` | GenomeManager, genome init/build flow, `get_metadata`, alias resolve/add/remove and alias trees | mixed |
| `managers/test_database.py` | `rgc.database.init`: folder creation, idempotency, empty config, no writes under home; `Refgenie` constructor coercion | unit |
| `managers/test_sequence.py` | `rgc.sequence.get` against a real RefgetStore, remote fallback for metadata-only genomes, `parse_locus` | unit |
| `managers/test_pull.py` | `rgc.transfer.pull` end to end, rollback, `pull_multiple` | mixed |
| `managers/test_transfer.py` | `TransferManager` (`rgc.transfer`) bulk skeleton with mocked dependencies: `pull_genomes` filtering and expansion, confirm refusal, `init_genomes`, `mirror` registration | unit |
| `managers/test_remote.py` | `RemoteManager` (`rgc.remote`): the id-or-name naming rule, add/upsert/remove, links, `unpushed`/`mark_pushed`, status, `pushed_urls`, push intent from `stage.create(push_to=...)` | component |
| `managers/test_servers.py` | `rgc.servers`: subscriptions, clients, remote listing, seekr, `resolve_alias`, `find_collection`; `genome.init_from_remote` | mixed |
| `managers/test_recipes.py` | Recipe CRUD and validation, multi-version asset classes and recipes, Snakefile generation | mixed |
| `managers/test_server_client.py` | `RefgenieserverClient`, store-mode selection, download progress | unit |
| `managers/test_serving_modes.py` | Serving-mode resolution and staging by mode, GA4GH DRS endpoints vs serving modes, remote download redirects (DRS and v4) | mixed |
| `managers/test_sources.py` | Genome sources: service-info bootstrap, source cache, rgstore detection, remote genomes, browse/sync | unit |
| `db/test_db.py` | ORM round-trips, Alembic chain, the remote-name migration up and down, manager root, database config and migration state, error hierarchy, fixtures-vs-helpers hygiene | mixed |
| `utils/test_build_provenance.py` | `build_digest`: the identity of a build, as opposed to its output | unit |
| `utils/test_dir_digest.py` | `get_dir_digest`: the asset primary key | unit |
| `utils/test_utils.py` | `refgenie.utils.*`: IO/YAML, checksums, encryption, tarballs, progress columns, confirmation prompts | mixed |
| `cli/test_cli.py` | argv parsing, model validation, help rendering, exit codes, dispatch, `refgenie id`, the `refgenie-build-fasta` console script, the model-field-is-consumed regression guard | mixed |
| `cli/test_remote_push.py` | The `refgenie remote status` handler, the push workflow, catalog export/import | mixed |
| `mcp/test_mcp.py` | MCP tools and the `mcp`-extras import guard | unit |
| `server/test_server.py` | The server HTTP app: endpoints, service info, data-channel router, route guards | mixed |
| `server/test_app.py` | `create_app(mode)` and mode isolation, the shared JSON surface, `/service-info` bootstrap, route/operationId uniqueness, the SPA-route/API non-collision guard, SPA serving, the OpenAPI-to-TypeScript / SKILL.md drift guards | unit |
| `server/test_bridge.py` | The localhost bridge: `/ping` contract, bridge modes, bridge CORS, the local-network-access preflight, cross-origin action policy | unit |
| `server/test_jobs.py` | JobManager lifecycle/cancel/coalescing/history, pull/build error mapping, the jobs HTTP API and its multiplexed SSE stream, log capture | unit |
| `server/test_local_actions.py` | The local actions API `/v1/actions/*`: envelope, security stack, build preflight, synchronous curation endpoints, route inventory | mixed |
| `server/test_real_jobs.py` | Real pulls and builds through the job machinery: threaded builds, stage validation first, one real pull through the manager, the assembled local app end to end | mixed |
| `test_publish_catalog.py` | Alias coverage in the published catalog export | component |

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
pytest -m "unit or component"             # unit + component
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
pytestmark = pytest.mark.component  # whole module
```
```python
@pytest.mark.component  # single class
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
against a test double (see `managers/test_genome.py::TestGenomeInitFailure`) stay `unit`.

## A note on `-s` and pypiper

`addopts` in `pyproject.toml` does not pass `-s`; the suite passes with
pytest's default capture. Add `-s` yourself for readable log output. It is not
required.

`PipelineManager.__init__` calls
`atexit.register(self._exit_handler)` and never unregisters. Every
asset-building test leaves one behind (measured: 37 queued after 30 tests), so
a full `component` run fires a few hundred of them at interpreter shutdown.
Harmless (~0.5s) but it explains the odd activity after the last test.

## Two homes: fixtures vs helpers

There is one rule about where shared test code lives, and the
fixture-hygiene tests in `tests/db/test_db.py` enforce it:

- **`tests/helpers.py`** — every plain importable function, class and constant.
  No fixtures. Its top-level imports are limited to the stdlib, `pytest` and
  `unittest.mock`; `refgenie`, `sqlmodel` and `fastapi` are imported inside the
  functions that need them, so the mock-only tests still run without the
  optional extras and `tests/e2e/cli_compat/` stays library-free.
- **`tests/conftest.py`** — fixtures and pytest hooks only. It imports from
  `tests/helpers.py` what its fixtures need. The same rule holds for the
  nested conftests: test modules import plain helpers from the `helpers.py`
  beside them (`tests/integration/helpers.py`, `tests/e2e/cli_compat/helpers.py`),
  never from a conftest.

Never redefine a fixture name the root conftest already owns. A local
`refgenie_built` once silently changed the asset name for ten tests, and a
duplicate `fixtures_path` shadowed the root one for the whole integration
suite; the fixture-hygiene tests in `db/test_db.py` now fail on either.

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
| `route_keys`, `effective_api_routes`, `assert_no_duplicate_routes`, `assert_unique_operation_ids` | FastAPI route inspection (the route-registration guards) |
| `run_cli`, `popen_cli`, `cli_argv`, `find_free_port`, `wait_for_server` | the one CLI subprocess launcher |
| `server_with_asset`, `seed_staged_asset`, `seed_asset`, `create_asset_class`, `create_asset_on_disk` | direct DB seeding for server tests |
| `make_server_client_world` | a server + client pair for pull/list-remote |
| `mock_server_client`, `real_client_no_init`, `MockRemoteSource`, `mocked_puller` | the pull/seekr mocks |
| `fasta_asset`, `only_asset`, `content_dir`, `build_dir`, `build_flag`, `boom` | the standard built world |
| `db_snapshot`, `disk_snapshot`, `assets_rows`, `asset_name_rows` | state snapshots |
| `add_asset_from_files`, `stage_copy`, `make_command_values` | asset construction |
| `TESTS_DATA_DIR`, `RCRSD_FASTA`, `DEMO_FASTA`, `GENOME`, `GROUP`, `ASSET`, `OMIT`, `PUBLIC_ORIGIN`, `EVIL_ORIGIN`, `CAPABILITY_KEY_SET` | constants |

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
stops testing the one that ships, so a route can be broken in the real app
without a single failing test.

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

The `Refgenie` class delegates to managers. Don't look for methods like `refgenie.get_aliases()` - use the manager attributes:

| Property | Manager | Common Methods |
|----------|---------|----------------|
| `refgenie.alias` | AliasManager | `resolve()`, `get_for_genome()`, `add()`, `remove()` |
| `refgenie.genome` | GenomeManager | `add()`, `get()`, `initialize_genome()` |
| `refgenie.asset` | AssetManager | `content.add()`, `get()`, `seek()`, `remove()` |
| `refgenie.servers` | ServerManager | `subscriptions()`, `subscribe()`, `client()`, `seek()`, `resolve_alias()`, `find_collection()` |
| `refgenie.sources` | SourceManager | `add_channel()`, `list_channels()` |
| `refgenie.remote` | RemoteManager | `add()`, `get()`, `link()`, `unpushed()`, `mark_pushed()`, `pushed_urls()` |

Example: `refgenie.alias.get_for_genome(genome_digest)` not `refgenie.get_aliases()`
