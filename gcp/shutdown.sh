#!/bin/bash
# GCP VM shutdown script for graceful container shutdown during preemption
#
# This script runs when the VM receives a preemption signal from GCP.
# GCP provides ~30 seconds for shutdown, so we use 25s timeout to stay safe.
#
# All output is logged to /var/log/shutdown-script.log for debugging.

set -euo pipefail

# Log all output to shutdown script log
exec > >(tee -a /var/log/shutdown-script.log)
exec 2>&1

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo_info() { echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"; }
echo_warn() { echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"; }
echo_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"; }

echo_info "Preemption detected, starting graceful shutdown..."

# Change to application directory
APP_DIR="/opt/qm-nmr-calc"
if [[ ! -d "$APP_DIR" ]]; then
    echo_warn "Application directory $APP_DIR not found, nothing to shut down"
    exit 0
fi

cd "$APP_DIR"

# Check if containers are running
if ! docker compose ps -q 2>/dev/null | grep -q .; then
    echo_info "No containers running, nothing to shut down"
    exit 0
fi

# Stop containers with 25-second timeout (stay within 30s preemption window)
# This overrides the worker's 300s stop_grace_period
echo_info "Stopping Docker Compose services (25s timeout)..."
if docker compose stop --timeout 25; then
    echo_info "Containers stopped gracefully"
else
    echo_warn "Some containers may not have stopped cleanly"
fi

echo_info "Shutdown script complete at $(date '+%Y-%m-%d %H:%M:%S')"
