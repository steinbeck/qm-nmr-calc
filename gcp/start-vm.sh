#!/bin/bash
# Start GCP Spot VM for qm-nmr-calc
#
# Starts a stopped VM. The VM retains its configuration,
# attached disks, and static IP address.
#
# Usage:
#   ./start-vm.sh

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

# Check if already running
if [[ "$CURRENT_STATUS" == "RUNNING" ]]; then
    echo_warn "VM '$VM_NAME' is already running"

    # Get IP address
    EXTERNAL_IP=$(gcloud compute instances describe "$VM_NAME" \
        --zone="$GCP_ZONE" \
        --format="value(networkInterfaces[0].accessConfigs[0].natIP)")

    echo ""
    echo "External IP: $EXTERNAL_IP"
    echo "Use ./status-vm.sh to see full VM status"
    exit 0
fi

# ============================================================================
# Start the VM
# ============================================================================
echo ""
echo_info "Starting VM '$VM_NAME'..."
echo_info "This will take 2-3 minutes while containers start..."
gcloud compute instances start "$VM_NAME" --zone="$GCP_ZONE" --quiet

# Wait a moment for the VM to get its IP
sleep 5

# Get the external IP
EXTERNAL_IP=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(networkInterfaces[0].accessConfigs[0].natIP)")

echo ""
echo "=============================================="
echo_info "VM Started Successfully!"
echo "=============================================="
echo ""
echo "VM Name:     $VM_NAME"
echo "Zone:        $GCP_ZONE"
echo "Status:      RUNNING"
echo "External IP: $EXTERNAL_IP"
echo ""
echo "Startup time:"
echo "  The VM needs 2-3 minutes for Docker containers to start."
echo "  Monitor progress with: ./logs-vm.sh"
echo ""
echo "Once ready, access your deployment at your configured domain."
echo ""
echo "Useful commands:"
echo "  ./status-vm.sh  - Check VM status"
echo "  ./ssh-vm.sh     - SSH into the VM"
echo "  ./logs-vm.sh    - View container logs"
echo ""
