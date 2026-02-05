#!/bin/bash
# Delete GCP Spot VM for qm-nmr-calc
#
# Deletes the VM but PRESERVES the persistent disk.
# This allows you to recreate the VM later without losing job data.
#
# To delete everything including data, use ./teardown-infrastructure.sh
#
# Usage:
#   ./delete-vm.sh

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
# Check VM exists
# ============================================================================
VM_NAME="${RESOURCE_PREFIX}-vm"
echo_info "Checking VM '$VM_NAME' in zone '$GCP_ZONE'..."

if ! gcloud compute instances describe "$VM_NAME" --zone="$GCP_ZONE" &>/dev/null; then
    echo_warn "VM '$VM_NAME' does not exist in zone '$GCP_ZONE'"
    echo "Nothing to delete."
    exit 0
fi

# Get current status for display
CURRENT_STATUS=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(status)")

echo_info "Found VM '$VM_NAME' with status: $CURRENT_STATUS"

# ============================================================================
# Confirmation prompt
# ============================================================================
echo ""
echo_warn "WARNING: This will delete VM '$VM_NAME'"
echo ""
echo "What will be DELETED:"
echo "  - VM instance (compute charges stop)"
echo "  - Boot disk (20GB)"
echo ""
echo "What will be PRESERVED:"
echo "  - Persistent disk '${RESOURCE_PREFIX}-data' (all job data)"
echo "  - Static IP '${RESOURCE_PREFIX}-ip'"
echo "  - Firewall rules"
echo ""
echo "You can recreate the VM later with ./deploy-vm.sh"
echo "All job data will be retained on the persistent disk."
echo ""

read -p "Are you sure you want to delete the VM? (yes/no): " CONFIRM
if [[ "$CONFIRM" != "yes" ]]; then
    echo_info "Delete cancelled"
    exit 0
fi

# ============================================================================
# Delete the VM
# ============================================================================
echo ""
echo_info "Deleting VM '$VM_NAME'..."
gcloud compute instances delete "$VM_NAME" --zone="$GCP_ZONE" --quiet

echo ""
echo "=============================================="
echo_info "VM Deleted Successfully!"
echo "=============================================="
echo ""
echo "Deleted:"
echo "  - VM instance: $VM_NAME"
echo "  - Boot disk: ${VM_NAME} (20GB)"
echo ""
echo "Preserved:"
echo "  - Persistent disk: ${RESOURCE_PREFIX}-data"
echo "  - Static IP: ${RESOURCE_PREFIX}-ip"
echo "  - Firewall rules"
echo ""
echo "To recreate the VM with your preserved data:"
echo "  ./deploy-vm.sh"
echo ""
echo "To delete everything including data:"
echo "  ./teardown-infrastructure.sh"
echo ""
