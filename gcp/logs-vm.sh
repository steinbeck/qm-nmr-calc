#!/bin/bash
# View container logs on GCP Spot VM for qm-nmr-calc
#
# Streams Docker Compose logs from the running containers.
#
# Usage:
#   ./logs-vm.sh              # Stream all container logs
#   ./logs-vm.sh api          # Stream only API container logs
#   ./logs-vm.sh worker       # Stream only worker container logs
#
# Press Ctrl+C to stop streaming logs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
echo_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
echo_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# ============================================================================
# Load configuration
# ============================================================================
source "$SCRIPT_DIR/lib/config.sh"
load_config "$SCRIPT_DIR/config.toml" || exit 1

# Validate project is set (load_config ensures this, but be explicit)
if [[ -z "${GCP_PROJECT_ID:-}" ]]; then
    echo_error "GCP_PROJECT_ID not set after loading config"
    exit 1
fi

# ============================================================================
# Check gcloud authentication
# ============================================================================
ACTIVE_ACCOUNT=$(gcloud auth list --filter=status:ACTIVE --format="value(account)" 2>/dev/null || true)
if [[ -z "$ACTIVE_ACCOUNT" ]]; then
    echo_error "Not authenticated with gcloud"
    echo "Run: gcloud auth login"
    echo "Then: gcloud config set project $GCP_PROJECT_ID"
    exit 1
fi

# Set project
gcloud config set project "$GCP_PROJECT_ID" --quiet

# ============================================================================
# Check VM exists and is running
# ============================================================================
VM_NAME="${RESOURCE_PREFIX}-vm"

# Detect VM zone dynamically (v2.7 selects zone at deployment time)
GCP_ZONE=$(gcloud compute instances list --project="$GCP_PROJECT_ID" \
    --filter="name=${VM_NAME}" --format="value(zone)" --quiet 2>/dev/null | head -1)
if [[ -z "$GCP_ZONE" ]]; then
    echo_error "VM '$VM_NAME' not found in any zone"
    echo "Create the VM first with ./deploy-auto.sh"
    exit 1
fi
echo_info "Found VM '$VM_NAME' in zone '$GCP_ZONE'"

# Get current status
CURRENT_STATUS=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(status)")

if [[ "$CURRENT_STATUS" != "RUNNING" ]]; then
    echo_error "VM '$VM_NAME' is not running (status: $CURRENT_STATUS)"
    echo ""
    echo "Start the VM first with: ./start-vm.sh"
    exit 1
fi

# ============================================================================
# Build logs command
# ============================================================================
SERVICE="${1:-}"
COMPOSE_CMD="docker compose -f /opt/qm-nmr-calc/docker-compose.yml -f /opt/qm-nmr-calc/docker-compose.gcp.yml logs -f"

if [[ -n "$SERVICE" ]]; then
    echo_info "Streaming logs for service '$SERVICE' on VM '$VM_NAME'..."
    COMPOSE_CMD="$COMPOSE_CMD $SERVICE"
else
    echo_info "Streaming all container logs on VM '$VM_NAME'..."
fi

echo_info "Press Ctrl+C to stop"
echo ""

# ============================================================================
# Stream logs
# ============================================================================
gcloud compute ssh "$VM_NAME" --zone="$GCP_ZONE" --command="$COMPOSE_CMD"
