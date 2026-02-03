#!/bin/bash
# Validation script for Docker Compose integration
# Tests stack startup, health checks, persistence, restart policies, and logs
#
# Usage:
#   ./scripts/validate-compose.sh           # Run all tests, cleanup at end
#   ./scripts/validate-compose.sh --no-cleanup  # Run tests, leave stack running
#   ./scripts/validate-compose.sh --help    # Show this help message
#
# Requirements:
#   - Docker and Docker Compose installed
#   - Run from project root directory

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
API_PORT="${API_PORT:-8000}"
HEALTH_TIMEOUT=120  # seconds to wait for health checks
RESTART_TIMEOUT=30  # seconds to wait for service restart

# Flags
CLEANUP=true

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --no-cleanup)
            CLEANUP=false
            shift
            ;;
        --help|-h)
            echo "Docker Compose Integration Validation"
            echo ""
            echo "Usage:"
            echo "  $0              Run all tests, cleanup at end"
            echo "  $0 --no-cleanup Run tests, leave stack running"
            echo "  $0 --help       Show this help message"
            echo ""
            echo "Tests performed:"
            echo "  1. Stack Startup      - Build and start services"
            echo "  2. API Health         - Verify API health endpoint"
            echo "  3. Worker Health      - Verify Huey consumer running"
            echo "  4. Volume Persistence - Data survives restart"
            echo "  5. Restart Policy     - Service auto-recovers from crash"
            echo "  6. Logs Accessible    - Can view service logs"
            echo ""
            echo "Environment variables:"
            echo "  API_PORT    Port to test (default: 8000)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Helper function to print test status
check_status() {
    local test_name="$1"
    local status="$2"

    if [ "$status" -eq 0 ]; then
        echo -e "${GREEN}[PASS]${NC} $test_name"
        return 0
    else
        echo -e "${RED}[FAIL]${NC} $test_name"
        return 1
    fi
}

# Helper function to print info messages
info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# Helper function to print warning messages
warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

# Cleanup function
cleanup() {
    if [ "$CLEANUP" = true ]; then
        info "Cleaning up: removing containers and volumes..."
        docker compose down -v 2>/dev/null || true
    else
        info "Leaving stack running (--no-cleanup specified)"
        info "To stop later: docker compose down -v"
    fi
}

# Set up cleanup on exit (only if CLEANUP is true)
trap 'if [ "$CLEANUP" = true ]; then cleanup; fi' EXIT

echo ""
echo "=========================================="
echo "Docker Compose Integration Validation"
echo "=========================================="
echo ""

# Track overall status
ALL_PASSED=true

# ============================================================================
# Test 1: Stack Startup
# ============================================================================
echo -e "\n${YELLOW}Test 1: Stack Startup${NC}"
info "Building and starting services (this may take a few minutes)..."

# Clean any existing containers first
docker compose down -v 2>/dev/null || true

# Build and start
if docker compose up -d --build 2>&1 | tail -5; then
    info "Services starting, waiting for health checks (timeout: ${HEALTH_TIMEOUT}s)..."
else
    check_status "Docker compose up" 1
    ALL_PASSED=false
fi

# Wait for services to be healthy
WAIT_COUNT=0
while [ $WAIT_COUNT -lt $HEALTH_TIMEOUT ]; do
    # Count healthy services - use tr to ensure single line output
    WORKER_HEALTHY=$(docker compose ps api worker --format json 2>/dev/null | grep -c '"Health":"healthy"' 2>/dev/null | tr -d '\n' || echo "0")
    # Ensure we have a valid number
    WORKER_HEALTHY=${WORKER_HEALTHY:-0}

    if [ "$WORKER_HEALTHY" -ge 2 ] 2>/dev/null; then
        break
    fi

    sleep 2
    WAIT_COUNT=$((WAIT_COUNT + 2))
    if [ $((WAIT_COUNT % 10)) -eq 0 ]; then
        info "Still waiting... (${WAIT_COUNT}s/${HEALTH_TIMEOUT}s)"
    fi
done

# Check final health status
if [ "$WORKER_HEALTHY" -ge 2 ] 2>/dev/null; then
    check_status "Both services healthy" 0
else
    check_status "Both services healthy (timeout after ${HEALTH_TIMEOUT}s)" 1
    ALL_PASSED=false
    docker compose ps
fi

# ============================================================================
# Test 2: API Health
# ============================================================================
echo -e "\n${YELLOW}Test 2: API Health${NC}"

HTTP_STATUS=$(curl -s -o /dev/null -w "%{http_code}" "http://localhost:${API_PORT}/health" 2>/dev/null || echo "000")

if [ "$HTTP_STATUS" = "200" ]; then
    check_status "API health endpoint returns 200" 0
else
    check_status "API health endpoint (got HTTP $HTTP_STATUS)" 1
    ALL_PASSED=false
fi

# ============================================================================
# Test 3: Worker Health
# ============================================================================
echo -e "\n${YELLOW}Test 3: Worker Health${NC}"

if docker compose exec -T worker pgrep -f huey_consumer > /dev/null 2>&1; then
    check_status "Huey consumer process running" 0
else
    check_status "Huey consumer process running" 1
    ALL_PASSED=false
fi

# ============================================================================
# Test 4: Volume Persistence
# ============================================================================
echo -e "\n${YELLOW}Test 4: Volume Persistence${NC}"

# Create a marker file
info "Creating marker file in data volume..."
docker compose exec -T api touch /app/data/test-marker 2>/dev/null

# Verify marker was created
if docker compose exec -T api ls /app/data/test-marker > /dev/null 2>&1; then
    info "Marker file created successfully"
else
    check_status "Create marker file" 1
    ALL_PASSED=false
fi

# Restart stack (down then up, without removing volumes)
info "Restarting stack (preserving volumes)..."
docker compose down 2>/dev/null
docker compose up -d 2>/dev/null

# Wait for services to be healthy again
info "Waiting for services to be healthy..."
WAIT_COUNT=0
while [ $WAIT_COUNT -lt $HEALTH_TIMEOUT ]; do
    WORKER_HEALTHY=$(docker compose ps api worker --format json 2>/dev/null | grep -c '"Health":"healthy"' 2>/dev/null | tr -d '\n' || echo "0")
    WORKER_HEALTHY=${WORKER_HEALTHY:-0}

    if [ "$WORKER_HEALTHY" -ge 2 ] 2>/dev/null; then
        break
    fi

    sleep 2
    WAIT_COUNT=$((WAIT_COUNT + 2))
done

# Verify marker still exists
if docker compose exec -T api ls /app/data/test-marker > /dev/null 2>&1; then
    check_status "Marker file persists after restart" 0
    # Clean up marker
    docker compose exec -T api rm /app/data/test-marker 2>/dev/null || true
else
    check_status "Marker file persists after restart" 1
    ALL_PASSED=false
fi

# ============================================================================
# Test 5: Restart Policy
# ============================================================================
echo -e "\n${YELLOW}Test 5: Restart Policy${NC}"

info "Simulating API crash (killing uvicorn)..."
docker compose exec -T api pkill -9 uvicorn 2>/dev/null || true

# Wait for restart
info "Waiting for service to auto-recover (timeout: ${RESTART_TIMEOUT}s)..."
sleep 5  # Give Docker time to detect the crash

WAIT_COUNT=0
RECOVERED=false
while [ $WAIT_COUNT -lt $RESTART_TIMEOUT ]; do
    HTTP_STATUS=$(curl -s -o /dev/null -w "%{http_code}" "http://localhost:${API_PORT}/health" 2>/dev/null || echo "000")

    if [ "$HTTP_STATUS" = "200" ]; then
        RECOVERED=true
        break
    fi

    sleep 2
    WAIT_COUNT=$((WAIT_COUNT + 2))
done

if [ "$RECOVERED" = true ]; then
    check_status "API auto-recovered after crash" 0
else
    check_status "API auto-recovered after crash (timeout after ${RESTART_TIMEOUT}s)" 1
    ALL_PASSED=false
fi

# ============================================================================
# Test 6: Logs Accessible
# ============================================================================
echo -e "\n${YELLOW}Test 6: Logs Accessible${NC}"

# Test API logs
API_LOGS=$(docker compose logs --tail 10 api 2>&1)
if [ -n "$API_LOGS" ] && [ ${#API_LOGS} -gt 20 ]; then
    check_status "API logs accessible (${#API_LOGS} chars)" 0
else
    check_status "API logs accessible" 1
    ALL_PASSED=false
fi

# Test worker logs
WORKER_LOGS=$(docker compose logs --tail 10 worker 2>&1)
if [ -n "$WORKER_LOGS" ] && [ ${#WORKER_LOGS} -gt 20 ]; then
    check_status "Worker logs accessible (${#WORKER_LOGS} chars)" 0
else
    check_status "Worker logs accessible" 1
    ALL_PASSED=false
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=========================================="
echo "Validation Summary"
echo "=========================================="

if [ "$ALL_PASSED" = true ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    echo ""
    docker compose ps
    exit 0
else
    echo -e "${RED}Some tests failed.${NC}"
    echo ""
    echo "Check logs with: docker compose logs"
    echo "Check status with: docker compose ps"
    exit 1
fi
