#!/bin/bash
# Test Services Management Script
# Manages services required for integration tests (PostgreSQL, HTTP data channel)
#
# Supports parallel test runs via unique container names and ports.
# Set REFGENIE_TEST_RUN_ID to share services across script invocations.
#
# Usage:
#   ./tests/scripts/services.sh start   # Start all services
#   ./tests/scripts/services.sh stop    # Stop all services
#   ./tests/scripts/services.sh status  # Show service status

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR/../.."
TEST_DATA_DIR="$PROJECT_ROOT/tests/data"

# Generate unique run ID if not provided (enables parallel test runs)
RUN_ID="${REFGENIE_TEST_RUN_ID:-$$}"

# Use environment variables or generate unique values
CONTAINER_NAME="${REFGENIE_TEST_CONTAINER:-refgenie-postgres-test-${RUN_ID}}"
DB_PORT="${REFGENIE_TEST_DB_PORT:-$(( 5433 + (RANDOM % 1000) ))}"
HTTP_PORT="${REFGENIE_TEST_HTTP_PORT:-$(( 18765 + (RANDOM % 1000) ))}"

DB_USER="testuser"
DB_PASS="testpass"
DB_NAME="refgenie_test"

# Export for child processes and test-integration.sh
export REFGENIE_TEST_RUN_ID="$RUN_ID"
export REFGENIE_TEST_CONTAINER="$CONTAINER_NAME"
export REFGENIE_TEST_DB_PORT="$DB_PORT"
export REFGENIE_TEST_HTTP_PORT="$HTTP_PORT"

# PID file for HTTP server
HTTP_PID_FILE="/tmp/refgenie-http-test-${RUN_ID}.pid"

start_postgres() {
    echo "Starting PostgreSQL..."
    echo "  Container: $CONTAINER_NAME"
    echo "  Port: $DB_PORT"

    # Remove existing container if it exists
    docker rm -f "$CONTAINER_NAME" 2>/dev/null || true

    # Start PostgreSQL container with tmpfs for speed
    docker run -d \
        --name "$CONTAINER_NAME" \
        -e POSTGRES_USER="$DB_USER" \
        -e POSTGRES_PASSWORD="$DB_PASS" \
        -e POSTGRES_DB="$DB_NAME" \
        -p "${DB_PORT}:5432" \
        --tmpfs /var/lib/postgresql/data \
        postgres:17

    # Wait for healthy status (up to 30 seconds)
    echo "Waiting for PostgreSQL to be ready..."
    for i in {1..30}; do
        if docker exec "$CONTAINER_NAME" pg_isready -U "$DB_USER" -d "$DB_NAME" 2>/dev/null; then
            echo "PostgreSQL is ready!"
            return 0
        fi
        sleep 1
    done

    echo "Failed to start PostgreSQL"
    docker logs "$CONTAINER_NAME"
    return 1
}

stop_postgres() {
    echo "Stopping PostgreSQL ($CONTAINER_NAME)..."
    docker rm -f "$CONTAINER_NAME" 2>/dev/null || true
}

start_http() {
    echo "Starting HTTP data channel server..."
    echo "  Port: $HTTP_PORT"
    echo "  Serving: $TEST_DATA_DIR"

    # Check if port is already in use
    if command -v lsof &>/dev/null && lsof -i ":$HTTP_PORT" &>/dev/null; then
        echo "Warning: Port $HTTP_PORT already in use, attempting to free it..."
        fuser -k "$HTTP_PORT/tcp" 2>/dev/null || true
        sleep 1
    fi

    # Start Python HTTP server
    cd "$TEST_DATA_DIR"
    python3 -m http.server "$HTTP_PORT" --bind localhost &>/dev/null &
    echo $! > "$HTTP_PID_FILE"
    cd "$PROJECT_ROOT"

    sleep 1
    if [ -f "$HTTP_PID_FILE" ] && kill -0 "$(cat "$HTTP_PID_FILE")" 2>/dev/null; then
        echo "HTTP server is ready! (PID: $(cat "$HTTP_PID_FILE"))"
        return 0
    else
        echo "Failed to start HTTP server"
        return 1
    fi
}

stop_http() {
    echo "Stopping HTTP data channel server..."
    if [ -f "$HTTP_PID_FILE" ]; then
        local pid=$(cat "$HTTP_PID_FILE")
        if kill -0 "$pid" 2>/dev/null; then
            kill "$pid" 2>/dev/null || true
            wait "$pid" 2>/dev/null || true
        fi
        rm -f "$HTTP_PID_FILE"
    fi
}

show_status() {
    echo "=== Test Services Status (Run ID: $RUN_ID) ==="
    echo ""
    echo "PostgreSQL:"
    if docker ps -f "name=$CONTAINER_NAME" --format "  Container: {{.Names}} | Port: {{.Ports}} | Status: {{.Status}}" | grep -q .; then
        docker ps -f "name=$CONTAINER_NAME" --format "  Container: {{.Names}} | Port: {{.Ports}} | Status: {{.Status}}"
    else
        echo "  Not running"
    fi
    echo ""
    echo "HTTP Data Channel:"
    if [ -f "$HTTP_PID_FILE" ] && kill -0 "$(cat "$HTTP_PID_FILE")" 2>/dev/null; then
        echo "  Running (PID: $(cat "$HTTP_PID_FILE"), Port: $HTTP_PORT)"
    else
        echo "  Not running"
    fi
}

print_exports() {
    # Print export commands for manual use
    echo ""
    echo "To use these services manually, run:"
    echo "  export REFGENIE_TEST_RUN_ID=\"$RUN_ID\""
    echo "  export REFGENIE_TEST_DB_PORT=\"$DB_PORT\""
    echo "  export REFGENIE_TEST_HTTP_PORT=\"$HTTP_PORT\""
    echo "  export TEST_DB_URL=\"postgresql+psycopg://${DB_USER}:${DB_PASS}@localhost:${DB_PORT}/${DB_NAME}\""
    echo "  export DATA_CHANNEL_URL=\"http://localhost:${HTTP_PORT}\""
}

case "$1" in
    start)
        echo "=== Starting Test Services (Run ID: $RUN_ID) ==="
        start_postgres
        start_http
        print_exports
        echo ""
        echo "=== All services started ==="
        ;;
    stop)
        echo "=== Stopping Test Services ==="
        stop_http
        stop_postgres
        echo "=== All services stopped ==="
        ;;
    restart)
        $0 stop
        $0 start
        ;;
    status)
        show_status
        ;;
    *)
        echo "Usage: $0 {start|stop|restart|status}"
        echo ""
        echo "Manages services required for integration tests."
        echo ""
        echo "Environment variables for parallel execution:"
        echo "  REFGENIE_TEST_RUN_ID    - Unique identifier (default: PID)"
        echo "  REFGENIE_TEST_DB_PORT   - PostgreSQL port (default: 5433 + random)"
        echo "  REFGENIE_TEST_HTTP_PORT - HTTP server port (default: 18765 + random)"
        echo "  REFGENIE_TEST_CONTAINER - Container name (default: refgenie-postgres-test-\$RUN_ID)"
        exit 1
        ;;
esac
