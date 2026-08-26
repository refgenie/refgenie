#!/bin/bash
# Run integration tests with ephemeral test services
#
# This script handles all setup and teardown automatically:
# 1. Starts test services (PostgreSQL, HTTP data channel)
# 2. Runs integration tests
# 3. Tears down all services on exit (even on failure)
#
# Usage: ./tests/scripts/test-integration.sh [pytest args]
# Example: ./tests/scripts/test-integration.sh -v -k "test_health"
#
# For manual service control (debugging):
#   ./tests/scripts/services.sh start
#   RUN_INTEGRATION_TESTS=true pytest tests/integration/
#   ./tests/scripts/services.sh stop

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR/../.."
SERVICES_SCRIPT="$SCRIPT_DIR/services.sh"

# Generate unique run ID for this test session (enables parallel test runs)
export REFGENIE_TEST_RUN_ID="$$"

# Let services.sh generate the ports, but we need them for cleanup
export REFGENIE_TEST_DB_PORT=$(( 5433 + (RANDOM % 1000) ))
export REFGENIE_TEST_HTTP_PORT=$(( 18765 + (RANDOM % 1000) ))
export REFGENIE_TEST_CONTAINER="refgenie-postgres-test-${REFGENIE_TEST_RUN_ID}"

SERVICES_STARTED=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Cleanup function - runs on exit (success or failure)
cleanup() {
    local exit_code=$?
    echo ""
    if [ "$SERVICES_STARTED" = true ]; then
        log_info "Cleaning up..."
        "$SERVICES_SCRIPT" stop
    fi

    if [ $exit_code -eq 0 ]; then
        log_info "Integration tests completed successfully!"
    else
        log_error "Integration tests failed with exit code: $exit_code"
    fi

    exit $exit_code
}

trap cleanup EXIT INT TERM

echo "=============================================="
echo "   Refgenie Integration Tests"
echo "   Run ID: $REFGENIE_TEST_RUN_ID"
echo "=============================================="
echo ""

# Check prerequisites
log_info "Checking prerequisites..."

if ! command -v docker &> /dev/null; then
    log_error "Docker is not installed or not in PATH"
    exit 1
fi

if ! docker info &> /dev/null; then
    log_error "Docker daemon is not running"
    exit 1
fi

# Start services
log_info "Starting test services..."
"$SERVICES_SCRIPT" start
SERVICES_STARTED=true

# Export environment variables for tests
export RUN_INTEGRATION_TESTS=true
export TEST_DB_URL="postgresql+psycopg://testuser:testpass@localhost:${REFGENIE_TEST_DB_PORT}/refgenie_test"
export DATA_CHANNEL_URL="http://localhost:${REFGENIE_TEST_HTTP_PORT}"

# Run integration tests
echo ""
log_info "Running integration tests..."
echo "----------------------------------------------"

cd "$PROJECT_ROOT"

# Activate bulker crate for bioinformatics tools (samtools, bowtie2, bwa)
BULKER_MANIFEST="$PROJECT_ROOT/bulker_manifest.yaml"
if command -v bulker &> /dev/null && [ -f "$BULKER_MANIFEST" ]; then
    log_info "Running tests inside bulker crate: $BULKER_MANIFEST"
    bulker exec "$BULKER_MANIFEST" -- python3 -m pytest tests/integration/ "$@"
else
    log_info "bulker not found; running tests without containerized tools (some may skip)"
    python3 -m pytest tests/integration/ "$@"
fi
TEST_EXIT_CODE=$?

echo "----------------------------------------------"
exit $TEST_EXIT_CODE
