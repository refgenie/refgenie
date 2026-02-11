# Cross-Package Integration Tests

Integration tests for the **refgenie + refgenconf + refgenieserver** stack. These exercise the full client-server chain: server endpoints, `listr()`, `pull()`, `seek()`, alias resolution, subscriptions, config persistence, and CLI commands.

## How to run

```bash
./tests/scripts/test-integration.sh           # Full run
./tests/scripts/test-integration.sh -v         # Verbose
./tests/scripts/test-integration.sh -k "test_pull"  # Filter by name
```

The runner script:
1. Verifies sibling repos exist (`../refgenconf`, `../refgenieserver`)
2. Installs all three packages as editable (`pip install -e`)
3. Installs test dependencies (fastapi, uvicorn, httpx)
4. Sets `RUN_INTEGRATION_TESTS=true` and runs pytest

Tests are **skipped by default** when running `pytest` from the repo root. They only run when `RUN_INTEGRATION_TESTS=true` is set, so they won't interfere with normal unit test runs.

## Prerequisites

All three sibling repos must be present at the same directory level:

```
some-parent/
├── refgenie/          # This repo
├── refgenconf/        # Required sibling
└── refgenieserver/    # Required sibling
```

## Architecture

The tests use a **real threaded uvicorn server** (not TestClient mocks) because `refgenconf.pull()` uses three different HTTP mechanisms (`requests.get`, `urlopen`, `urlretrieve`) that would each need separate mocking.

**conftest.py** creates:
- A minimal FastAPI app mimicking refgenieserver's key endpoints with matching operation IDs
- A threaded uvicorn server on a random port
- Client `RefGenConf` instances subscribed to the test server
- Archived genome assets (`.tgz` files with checksums)

**test_server_client.py** has 6 test groups (33 tests total):

| Group | Tests | Coverage |
|-------|-------|----------|
| Server Endpoints | 11 | Direct HTTP verification of all endpoints |
| Client Operations | 6 | `listr()`, `pull()`, `seek()` with seek keys |
| Alias Operations | 4 | Digest resolution, custom aliases, aliases property |
| Subscription | 3 | Subscribe, unsubscribe, deduplication |
| Config Persistence | 2 | Write/reload cycle, pull persists to disk |
| CLI Integration | 4 | `refgenie list`, `listr`, `pull`, `seek` via subprocess |

## Key gotcha: GENOME_DIGEST

`GENOME_DIGEST` in `conftest.py` **must** match the digest that `SeqColClient.load_fasta()` computes from the test `FASTA_CONTENT`. During `pull()`, `initialize_genome()` recomputes the digest from the downloaded FASTA and calls `set_genome_alias(overwrite=True)`. If the server-reported digest doesn't match the computed one, two genome entries are created and alias resolution breaks with `UndefinedAliasError`.

If you change `FASTA_CONTENT`, recompute the digest:

```python
from refgenconf.seqcol import SeqColClient
import tempfile, os

with tempfile.NamedTemporaryFile(suffix='.fa', mode='w', delete=False) as f:
    f.write(FASTA_CONTENT)
    path = f.name

ssc = SeqColClient({})
digest, _ = ssc.load_fasta(path, gzipped=False)
print(digest)  # Use this as GENOME_DIGEST
os.unlink(path)
```
