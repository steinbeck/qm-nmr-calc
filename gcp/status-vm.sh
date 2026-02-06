#!/bin/bash
# Show GCP Spot VM status for qm-nmr-calc
#
# Displays VM status, IP address, machine type, and access information.
#
# Usage:
#   ./status-vm.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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
    echo ""
    echo "=============================================="
    echo_warn "VM Not Found"
    echo "=============================================="
    echo ""
    echo "VM Name:  $VM_NAME"
    echo "Status:   NOT CREATED"
    echo ""
    echo "To create a VM, run:"
    echo "  ./deploy-auto.sh"
    echo ""
    exit 0
fi
echo_info "Found VM '$VM_NAME' in zone '$GCP_ZONE'"

# ============================================================================
# Get VM details
# ============================================================================
VM_STATUS=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(status)")

MACHINE_TYPE=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(machineType)" | rev | cut -d'/' -f1 | rev)

EXTERNAL_IP=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(networkInterfaces[0].accessConfigs[0].natIP)" 2>/dev/null || echo "")

CREATION_TIME=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(creationTimestamp)")

# ============================================================================
# Display status
# ============================================================================
echo ""
echo "=============================================="
echo "  VM Status"
echo "=============================================="
echo ""
echo "VM Name:       $VM_NAME"
echo "Zone:          $GCP_ZONE"
echo "Machine Type:  $MACHINE_TYPE"
echo "Created:       $CREATION_TIME"
echo ""

# Status with color
case "$VM_STATUS" in
    RUNNING)
        echo -e "Status:        ${GREEN}$VM_STATUS${NC}"
        ;;
    TERMINATED|STOPPED)
        echo -e "Status:        ${YELLOW}$VM_STATUS${NC}"
        ;;
    STAGING|PROVISIONING|STOPPING)
        echo -e "Status:        ${BLUE}$VM_STATUS${NC}"
        ;;
    *)
        echo -e "Status:        ${RED}$VM_STATUS${NC}"
        ;;
esac

# External IP
if [[ -n "$EXTERNAL_IP" ]]; then
    echo "External IP:   $EXTERNAL_IP"
else
    echo "External IP:   N/A (VM is stopped)"
fi

echo ""

# ============================================================================
# Access information (if running)
# ============================================================================
if [[ "$VM_STATUS" == "RUNNING" ]]; then
    echo "=============================================="
    echo "  Access Information"
    echo "=============================================="
    echo ""
    if [[ -n "$EXTERNAL_IP" ]]; then
        echo "Web UI:  http://$EXTERNAL_IP"
        echo ""
    fi
    echo "SSH:"
    echo "  ./ssh-vm.sh"
    echo "  # or: gcloud compute ssh $VM_NAME --zone=$GCP_ZONE"
    echo ""
    echo "Container logs:"
    echo "  ./logs-vm.sh"
    echo ""
    echo "Quick commands:"
    echo "  ./stop-vm.sh   - Stop VM (pause billing)"
    echo "  ./ssh-vm.sh    - SSH into VM"
    echo "  ./logs-vm.sh   - View container logs"
    echo ""
elif [[ "$VM_STATUS" == "TERMINATED" || "$VM_STATUS" == "STOPPED" ]]; then
    echo "=============================================="
    echo "  VM is Stopped"
    echo "=============================================="
    echo ""
    echo "The VM is not running. No compute charges are being incurred."
    echo ""
    echo "To start the VM:"
    echo "  ./start-vm.sh"
    echo ""
    echo "To delete the VM (preserves data):"
    echo "  ./delete-vm.sh"
    echo ""
fi
