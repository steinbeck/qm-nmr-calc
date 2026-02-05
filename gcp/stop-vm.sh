#!/bin/bash
# Stop GCP Spot VM for qm-nmr-calc
#
# Stops the running VM to reduce costs when not in use.
# The persistent disk and static IP are preserved.
#
# Usage:
#   ./stop-vm.sh

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
if [[ ! -f "config.sh" ]]; then
    echo_error "config.sh not found!"
    echo "Copy config.sh.example to config.sh and set your GCP_PROJECT_ID:"
    echo "  cp config.sh.example config.sh"
    echo "  # Edit config.sh with your values"
    exit 1
fi

source ./config.sh

# Validate required variables
if [[ -z "${GCP_PROJECT_ID:-}" || "$GCP_PROJECT_ID" == "your-project-id" ]]; then
    echo_error "GCP_PROJECT_ID is not set or still has default value"
    echo "Edit config.sh and set your actual GCP project ID"
    exit 1
fi

# ============================================================================
# Check gcloud authentication
# ============================================================================
echo_info "Checking gcloud authentication..."
ACTIVE_ACCOUNT=$(gcloud auth list --filter=status:ACTIVE --format="value(account)" 2>/dev/null || true)
if [[ -z "$ACTIVE_ACCOUNT" ]]; then
    echo_error "Not authenticated with gcloud"
    echo "Run: gcloud auth login"
    echo "Then: gcloud config set project $GCP_PROJECT_ID"
    exit 1
fi
echo_info "Authenticated as: $ACTIVE_ACCOUNT"

# Set project
gcloud config set project "$GCP_PROJECT_ID" --quiet

# ============================================================================
# Check VM exists and get status
# ============================================================================
VM_NAME="${RESOURCE_PREFIX}-vm"
echo_info "Checking VM '$VM_NAME' in zone '$GCP_ZONE'..."

if ! gcloud compute instances describe "$VM_NAME" --zone="$GCP_ZONE" &>/dev/null; then
    echo_error "VM '$VM_NAME' not found in zone '$GCP_ZONE'"
    echo "Create the VM first with ./deploy-vm.sh"
    exit 1
fi

# Get current status
CURRENT_STATUS=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(status)")

echo_info "Current VM status: $CURRENT_STATUS"

# Check if already stopped
if [[ "$CURRENT_STATUS" == "TERMINATED" ]]; then
    echo_warn "VM '$VM_NAME' is already stopped"
    echo "No action needed. Use ./start-vm.sh to start the VM."
    exit 0
fi

# ============================================================================
# Stop the VM
# ============================================================================
echo ""
echo_info "Stopping VM '$VM_NAME'..."
gcloud compute instances stop "$VM_NAME" --zone="$GCP_ZONE" --quiet

echo ""
echo "=============================================="
echo_info "VM Stopped Successfully!"
echo "=============================================="
echo ""
echo "VM Name:    $VM_NAME"
echo "Zone:       $GCP_ZONE"
echo "Status:     TERMINATED (stopped)"
echo ""
echo "Billing note:"
echo "  - You are NO longer billed for compute time"
echo "  - Storage costs continue for boot disk and persistent disk"
echo "  - Static IP continues to incur charges"
echo ""
echo "To start the VM again:"
echo "  ./start-vm.sh"
echo ""
