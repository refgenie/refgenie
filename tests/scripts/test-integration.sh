#!/bin/bash
# Cross-package integration test runner for refgenie + refgenconf + refgenieserver.
#
# This script installs all three sibling packages from on-disk repos (not PyPI)
# and runs the integration tests using the REAL refgenieserver via app_factory.
#
# Usage:
#   ./tests/scripts/test-integration.sh [pytest-args...]
#
# Example:
#   ./tests/scripts/test-integration.sh -v
#   ./tests/scripts/test-integration.sh -k "test_pull"

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR/../.."
WORKSPACE="$(cd "$REPO_ROOT/.." && pwd)"

echo "=== Cross-package integration test runner ==="
echo "Repo root: $REPO_ROOT"
echo "Workspace: $WORKSPACE"

# Verify sibling repos exist
for pkg in refgenconf refgenieserver; do
    if [ ! -d "$WORKSPACE/$pkg" ]; then
        echo "ERROR: Sibling repo '$WORKSPACE/$pkg' not found."
        echo "Integration tests require all three packages to be present as sibling repos."
        exit 1
    fi
done

echo ""
echo "=== Installing sibling packages from on-disk repos ==="
pip install -e "$WORKSPACE/refgenconf" --quiet
pip install -e "$WORKSPACE/refgenieserver" --quiet
pip install -e "$REPO_ROOT" --quiet

# Install test dependencies
pip install uvicorn httpx fastapi pyyaml --quiet

echo ""
echo "=== Verifying imports ==="
python3 -c "from refgenieserver.app_factory import create_app; print('app_factory OK')"

echo ""
echo "=== Running integration tests ==="
export RUN_INTEGRATION_TESTS=true
cd "$REPO_ROOT"

# Run test files separately to avoid module-level global cross-contamination
# (both test files use create_app() which mutates refgenieserver.main globals)
echo "--- Server client integration tests ---"
python3 -m pytest tests/integration/test_server_client.py "$@"

echo ""
echo "--- Build integration tests ---"
python3 -m pytest tests/integration/test_build_integration.py "$@"

echo ""
echo "--- Build map/reduce mode tests ---"
python3 -m pytest tests/integration/test_build_map_mode.py "$@"
