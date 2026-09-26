## Contributing

We welcome contributions from the community.

## Development setup

1. Clone the repository:

```bash
git clone https://github.com/refgenie/refgenie1.git
cd refgenie1
```

2. Install in development mode:

```bash
pip install -e ".[dev,server,dash]"
```

3. Verify the installation:

```bash
refgenie --help
```

## Running tests

### Unit tests

Unit tests use in-memory SQLite and require no external dependencies. They run in about 7 seconds:

```bash
pytest tests/ --ignore=tests/integration
```

### Integration tests

Integration tests require Docker (for PostgreSQL) and take about 110 seconds. Use the provided script, which manages the PostgreSQL container lifecycle, HTTP data channel, and bulker crate activation:

```bash
./tests/scripts/test-integration.sh
```

Do not run `pytest tests/integration/` directly -- all tests will skip because the required environment variable and services will not be set up.

For debugging individual integration tests:

```bash
./tests/scripts/services.sh start
RUN_INTEGRATION_TESTS=true pytest tests/integration/ -k "test_name"
./tests/scripts/services.sh stop
```

## New assets

We are interested in working with groups who want to add assets into the refgenie system. Refgenie now uses a flexible asset class and recipe system -- you can define custom asset types and build recipes via YAML files without modifying refgenie's source code. See [Flexible asset types](flexible_asset_types.md) and [Data channels](data_channels.md) for details.

## Project structure

| Directory | Contents |
|---|---|
| `refgenie/` | Core library code |
| `refgenie/cli/` | CLI command definitions, handlers, and argument models |
| `refgenie/server/` | Built-in server (FastAPI routers, templates) |
| `refgenie/dash/` | Dashboard web UI |
| `refgenie/db/` | Database models (SQLModel) |
| `refgenie/managers/` | Asset, configuration, and genome managers |
| `refgenie/config/` | Configuration and database backend setup |
| `tests/` | Unit tests |
| `tests/integration/` | Integration tests (require Docker) |

## Suggestions and feedback

Please open an issue on the [GitHub issue tracker](https://github.com/refgenie/refgenie1/issues) with suggestions, bug reports, or other feedback.
