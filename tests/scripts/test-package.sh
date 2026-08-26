#!/bin/bash
# Build and validate the refgenie package in isolation -- mirrors CI
# (.github/workflows/package.yml) exactly, so a packaging failure can be
# reproduced locally without pushing.
#
# This script catches packaging bugs that unit tests can't detect:
# - A stale or missing web UI bundle riding into (or missing from) the wheel
# - Missing package data in general
# - Import errors in installed (non-editable) mode
#
# Usage: ./tests/scripts/test-package.sh
# Runtime: ~30-60 seconds (longer on the first run, while npm populates its cache)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR/../.."
cd "$PROJECT_ROOT"

# The app-factory smoke test below calls create_app(mode="local") with no
# refgenie_instance, so it resolves the caller's REAL refgenie config -- same
# as `refgenie dash` would. Point it at a disposable temp dir instead of
# ~/.refgenie: without this, a maintainer with $REFGENIE / $REFGENIE_DB_CONFIG_PATH
# set in their shell (or just a real ~/.refgenie) has this script open their
# actual catalog, which can fail on an unrelated migration mismatch and, worse,
# would let a bug in the smoke test itself write into real data.
TEMP_REFGENIE_HOME=$(mktemp -d)
export REFGENIE_HOME_PATH="$TEMP_REFGENIE_HOME"
unset REFGENIE REFGENIE_DB_CONFIG_PATH REFGENIE_GENOME_FOLDER REFGENIE_GENOME_STAGE_FOLDER REFGENIE_REFGET_STORE_URL

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_test() {
    echo -e "  ${GREEN}[TEST]${NC} $1"
}

# Cleanup function - runs on exit (success or failure)
cleanup() {
    local exit_code=$?
    echo ""
    log_info "Cleaning up..."

    if [ -d "$TEMP_VENV" ]; then
        rm -rf "$TEMP_VENV"
    fi

    if [ -n "$TEMP_REFGENIE_HOME" ] && [ -d "$TEMP_REFGENIE_HOME" ]; then
        rm -rf "$TEMP_REFGENIE_HOME"
    fi

    if [ $exit_code -eq 0 ]; then
        log_info "Package validation completed successfully!"
    else
        log_error "Package validation failed with exit code: $exit_code"
    fi

    exit $exit_code
}

trap cleanup EXIT INT TERM

echo "=============================================="
echo "   Refgenie Package Validation"
echo "=============================================="
echo ""

log_info "Checking prerequisites..."

if ! command -v python3 &> /dev/null; then
    log_error "python3 is not installed or not in PATH"
    exit 1
fi

if ! command -v uv &> /dev/null; then
    log_error "uv is not installed or not in PATH (https://docs.astral.sh/uv/)"
    exit 1
fi

# --- Step 1: the frontend bundle -------------------------------------------
#
# uv_build ships whatever is sitting in refgenie/server/webui/ at build time.
# A missing bundle must not silently produce a UI-less wheel: build it now if
# npm is available, or fail with the exact fix otherwise.
if [ ! -f "refgenie/server/webui/index.html" ]; then
    if command -v npm &> /dev/null; then
        log_info "refgenie/server/webui/index.html not found; building the frontend..."
        (cd frontend && npm ci && npm run build)
    else
        log_error "refgenie/server/webui/index.html is missing and npm is not on PATH."
        log_error "Install Node (see frontend/.node-version) and run:"
        log_error "    cd frontend && npm ci && npm run build"
        log_error "See docs/development.md for the full web UI dev workflow."
        exit 1
    fi
else
    log_info "Web UI bundle already present at refgenie/server/webui/index.html"
fi

# --- Step 2: build ------------------------------------------------------
log_info "Building wheel and sdist with 'uv build'..."
rm -rf dist
uv build

WHEEL_FILE=$(ls dist/*.whl 2>/dev/null | head -1)
SDIST_FILE=$(ls dist/*.tar.gz 2>/dev/null | head -1)
if [ -z "$WHEEL_FILE" ] || [ -z "$SDIST_FILE" ]; then
    log_error "uv build did not produce both a wheel and an sdist in dist/"
    exit 1
fi
log_info "Built: $(basename "$WHEEL_FILE"), $(basename "$SDIST_FILE")"

# --- Step 3: the wheel-content guard ----------------------------------------
#
# No --require-fresh here: that check is for publish.yml, where $GITHUB_SHA
# and a clean checkout are guaranteed. A local build is expected to be dirty.
log_info "Checking wheel and sdist carry the web UI bundle..."
python3 tests/scripts/check-wheel-assets.py "$WHEEL_FILE" "$SDIST_FILE"

# --- Step 4: install into a throwaway venv ----------------------------------
TEMP_VENV=$(mktemp -d)
log_info "Installing the wheel (non-editable, [dash] extras) into $TEMP_VENV..."
uv venv "$TEMP_VENV" --quiet
uv pip install --python "$TEMP_VENV/bin/python" "${WHEEL_FILE}[dash]" --quiet

# --- Step 5: smoke test against the INSTALLED package -----------------------
#
# This is the step that catches package-data bugs an editable checkout hides:
# it only ever sees refgenie/server/webui/ through the wheel's own contents.
echo ""
log_info "Running smoke tests against the installed package..."
echo "----------------------------------------------"

log_test "CLI module import..."
"$TEMP_VENV/bin/python" -c "from refgenie.cli.main import main; print('  OK: CLI module imports successfully')"

log_test "App factory serves the web UI (installed wheel)..."
"$TEMP_VENV/bin/python" -c "
from fastapi.testclient import TestClient
from refgenie.server.main import create_app, create_local_app
from refgenie.server.spa import resolve_web_dist, read_build_info

assert callable(create_local_app)
dist = resolve_web_dist()
assert dist is not None, 'web UI bundle missing from installed wheel'

app = create_app(mode='local')
# base_url: the local app's Host-header guard admits loopback names only,
# so the default 'testserver' host would 421 everything (see tests/test_web_ui.py).
with TestClient(app, base_url='http://localhost') as c:
    response = c.get('/')
    assert response.status_code == 200, response.status_code
    assert 'text/html' in response.headers['content-type']

print('  OK: create_app(mode=\"local\") builds and serves / from the installed wheel')
print(f'  Build info: {read_build_info(dist)}')
"

log_test "CLI --help..."
"$TEMP_VENV/bin/refgenie" --help > /dev/null 2>&1
echo "  OK: refgenie --help runs successfully"

log_test "Core package imports..."
"$TEMP_VENV/bin/python" -c "
from refgenie import Refgenie
from refgenie.models import BuildParams
from refgenie.db.tables import Asset, Genome
print('  OK: Core modules import successfully')
"

log_test "Package version..."
VERSION=$("$TEMP_VENV/bin/python" -c "from importlib.metadata import version; print(version('refgenie'))")
echo "  OK: Version $VERSION"

echo "----------------------------------------------"
log_info "All smoke tests passed!"

exit 0
