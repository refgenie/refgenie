# Refgenie Package

## Conventions

Repeated module basenames are fine and deliberate (`manager.py`, `const.py`,
`queries.py`, `helpers.py`), scoped by their package. One exception:
**`refgenie/models.py` is the only `models.py`** — it holds the domain types.
Wire models (request bodies, response envelopes, SSE frames) go in a
`schemas.py`; internal records that never hit the wire stay in the module that
owns them. `tests/test_layering.py` enforces this, and `docs/design-notes.md`
explains why. That file is the place to look for, and record, decisions that
are not obvious from the code.

## Testing

The suite is split into four tiers. Running everything takes about three
minutes; the inner loop takes about half a minute.

| Tier | Command | Needs |
|------|---------|-------|
| `unit` | `pytest` | nothing |
| `component` | `pytest -m component` | nothing |
| `e2e` | `pytest -m e2e` | `samtools` on PATH (some tests skip without it) |
| `integration` | `./tests/scripts/test-integration.sh` | Docker, bulker |

Test counts and wall times are not recorded here — they go stale. Run
`pytest --collect-only -q | tail -1` for a count and `pytest --durations=10`
to see where the time goes.

`pytest` with no arguments runs **unit only** (plus the integration tests, which
self-skip). `pytest -m "unit or component"` reproduces the old default scope. To
run literally everything except integration: `pytest -m "unit or component or e2e"`.

Tiers are assigned automatically by location (see `pytest_collection_modifyitems`
in `tests/conftest.py`) and can be overridden per module or per class with
`pytestmark = pytest.mark.component`. The dividing line: a test is `component`,
not `unit`, if it builds real genome folders, asset files or `.tgz` archives on
disk. See `tests/README.md` for which tier to write in.

IMPORTANT: "Run the tests" or "run integration tests" means `./tests/scripts/test-integration.sh`. Do NOT run `pytest tests/integration/` directly — all tests will skip because the env var and services won't be set up. The script handles PostgreSQL container lifecycle, HTTP data channel, bulker crate activation, and cleanup. There should be zero skipped tests when run this way.

**Manual service control (for debugging):**
```bash
./tests/scripts/services.sh start   # Start services
RUN_INTEGRATION_TESTS=true pytest tests/integration/ -k "test_name"
./tests/scripts/services.sh stop    # Stop services
```

See `tests/README.md` for detailed test organization and patterns.

### Server Endpoint Testing

Build the app with `create_app(refgenie_instance=test_rgc)` and use it directly.
`create_app(mode="local", ...)` builds the `refgenie dash` app from the same
factory; `make_server_app` / `make_local_app` in `tests/helpers.py` wrap both.
Never hand-build a `FastAPI()` from routers to get a TestClient.
`create_app` installs `app.dependency_overrides[get_refgenie] = lambda: rg`, and
`get_db_session` takes `get_refgenie` through `Depends`, so a handler's `rgc` and
its `session` are guaranteed to be the same database:

```python
from refgenie.server.main import create_app

app = create_app(refgenie_instance=test_rgc)
with TestClient(app, raise_server_exceptions=True) as client:
    response = client.get("/v4/genomes")
```

Do **not** patch `refgenie.server.dependencies._refgenie_instance` or
`server_main.refgenie`. That workaround existed only because `get_db_session`
used to call `get_refgenie()` directly, which FastAPI's override mechanism
cannot intercept. It no longer does.

Every server router — including GA4GH DRS — takes its dependencies through
`Depends`; there are no module-level globals to patch. HTML comes from the React
SPA in `frontend/`, built into `refgenie/server/webui/` and served by
`refgenie/server/spa.py`; there are no Jinja templates in the web layer.

To mock a **remote refgenieserver that a client pulls from**, wrap a second
Refgenie instance in `create_app` and inject the TestClient into the client's
source manager via `RefgenieserverClient(url, http_client=tc)` — helpers
`make_server_app(rg)` and `serve_refgenie(client_rg, server_rg, *urls)` live in
`tests/helpers.py` (plain functions live there; `tests/conftest.py` holds only
fixtures and hooks). Hand-rolled route fakes are not acceptable; the single
documented exception is `TestDownloadWithProgress` in `tests/test_server_client.py`
(the real FileResponse archive route cannot produce a no-Content-Length
response, which that client progress-bar regression test requires).
