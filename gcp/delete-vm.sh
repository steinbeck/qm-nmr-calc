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

# Detect VM zone dynamically (v2.7 selects zone at deployment time)
GCP_ZONE=$(gcloud compute instances list --project="$GCP_PROJECT_ID" \
    --filter="name=${VM_NAME}" --format="value(zone)" --quiet 2>/dev/null | head -1)
if [[ -z "$GCP_ZONE" ]]; then
    echo_warn "VM '$VM_NAME' does not exist in any zone"
    echo "Nothing to delete."
    exit 0
fi
echo_info "Found VM '$VM_NAME' in zone '$GCP_ZONE'"

# Get current status for display
CURRENT_STATUS=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(status)")

echo_info "Current VM status: $CURRENT_STATUS"

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
echo "You can recreate the VM later with ./deploy-auto.sh"
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
echo "  ./deploy-auto.sh"
echo ""
echo "To delete everything including data:"
echo "  ./teardown-infrastructure.sh"
echo ""
