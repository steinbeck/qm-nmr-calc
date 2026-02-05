#!/bin/bash
# SSH into GCP Spot VM for qm-nmr-calc
#
# Opens an SSH session to the VM, or executes a command.
#
# Usage:
#   ./ssh-vm.sh                    # Interactive SSH session
#   ./ssh-vm.sh "ls -la"           # Execute command and exit
#   ./ssh-vm.sh "docker ps"        # Check running containers

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

if ! gcloud compute instances describe "$VM_NAME" --zone="$GCP_ZONE" &>/dev/null; then
    echo_error "VM '$VM_NAME' not found in zone '$GCP_ZONE'"
    echo "Create the VM first with ./deploy-vm.sh"
    exit 1
fi

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
# SSH into VM
# ============================================================================
if [[ $# -eq 0 ]]; then
    # Interactive SSH session
    echo_info "Connecting to VM '$VM_NAME'..."
    gcloud compute ssh "$VM_NAME" --zone="$GCP_ZONE"
else
    # Execute command
    echo_info "Running command on VM '$VM_NAME'..."
    gcloud compute ssh "$VM_NAME" --zone="$GCP_ZONE" --command="$*"
fi
