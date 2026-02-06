#!/bin/bash
# Teardown GCP infrastructure for qm-nmr-calc
#
# Removes:
#   - Persistent disk (WARNING: deletes all job data!)
#   - Firewall rules
#   - Static external IP
#
# Usage:
#   ./teardown-infrastructure.sh          # Interactive with confirmation
#   ./teardown-infrastructure.sh --force  # Skip confirmation (use with caution!)

set -uo pipefail

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

# Parse arguments
FORCE=false
for arg in "$@"; do
    case $arg in
        --force|-f)
            FORCE=true
            shift
            ;;
    esac
done

# ============================================================================
# Load configuration
# ============================================================================
source "$SCRIPT_DIR/lib/config.sh"
load_config "$SCRIPT_DIR/config.toml" || exit 1

if [[ -z "${GCP_PROJECT_ID:-}" ]]; then
    echo_error "GCP_PROJECT_ID not set after loading config"
    exit 1
fi

# Set project
gcloud config set project "$GCP_PROJECT_ID" --quiet

# ============================================================================
# Detect zone and region
# ============================================================================
VM_NAME="${RESOURCE_PREFIX}-vm"
DISK_NAME="${RESOURCE_PREFIX}-data"

# Try to detect zone from existing VM first, then from disk
GCP_ZONE=$(gcloud compute instances list --project="$GCP_PROJECT_ID" \
    --filter="name=${VM_NAME}" --format="value(zone)" --quiet 2>/dev/null | head -1)

if [[ -z "$GCP_ZONE" ]]; then
    # VM might already be deleted, try to find zone from the persistent disk
    GCP_ZONE=$(gcloud compute disks list --project="$GCP_PROJECT_ID" \
        --filter="name=${DISK_NAME}" --format="value(zone)" --quiet 2>/dev/null | head -1)
fi

if [[ -z "$GCP_ZONE" ]]; then
    echo_warn "Could not detect zone from VM or disk (resources may already be deleted)"
    echo_warn "Proceeding with firewall rules and IP cleanup only..."
    ZONE_AVAILABLE=false
else
    # Derive region from zone (e.g., us-central1-a -> us-central1)
    GCP_REGION="${GCP_ZONE%-*}"
    ZONE_AVAILABLE=true
    echo_info "Detected zone: $GCP_ZONE"
    echo_info "Derived region: $GCP_REGION"
fi

echo_info "Teardown configuration:"
echo "  Project:  $GCP_PROJECT_ID"
if [[ "$ZONE_AVAILABLE" == "true" ]]; then
    echo "  Region:   $GCP_REGION"
    echo "  Zone:     $GCP_ZONE"
fi
echo "  Prefix:   $RESOURCE_PREFIX"
echo ""

# Confirmation prompt
if [[ "$FORCE" != "true" ]]; then
    echo ""
    echo_warn "WARNING: This will delete the following resources:"
    echo "  - Persistent disk '${RESOURCE_PREFIX}-data' with ALL job data"
    echo "  - Firewall rules for HTTP and SSH"
    echo "  - Static IP '${RESOURCE_PREFIX}-ip'"
    echo ""
    read -p "Are you sure you want to proceed? (yes/no): " CONFIRM
    if [[ "$CONFIRM" != "yes" ]]; then
        echo_info "Teardown cancelled"
        exit 0
    fi
fi

echo ""
echo_info "Starting infrastructure teardown..."

# ============================================================================
# Delete persistent disk (requires zone)
# ============================================================================
if [[ "$ZONE_AVAILABLE" == "true" ]]; then
    DISK_NAME="${RESOURCE_PREFIX}-data"
    echo ""
    echo_info "Deleting persistent disk '$DISK_NAME'..."
    if gcloud compute disks describe "$DISK_NAME" --zone="$GCP_ZONE" &>/dev/null; then
        gcloud compute disks delete "$DISK_NAME" --zone="$GCP_ZONE" --quiet || echo_warn "Failed to delete disk (may be attached to VM)"
        echo_info "Persistent disk '$DISK_NAME' deleted"
    else
        echo_warn "Persistent disk '$DISK_NAME' does not exist, skipping"
    fi
else
    echo ""
    echo_warn "Skipping disk deletion (zone unknown)"
fi

# ============================================================================
# Delete firewall rules
# ============================================================================
echo ""
echo_info "Deleting firewall rules..."

HTTP_RULE="${RESOURCE_PREFIX}-allow-http"
SSH_RULE="${RESOURCE_PREFIX}-allow-ssh"

for RULE in "$HTTP_RULE" "$SSH_RULE"; do
    if gcloud compute firewall-rules describe "$RULE" &>/dev/null; then
        gcloud compute firewall-rules delete "$RULE" --quiet || echo_warn "Failed to delete firewall rule '$RULE'"
        echo_info "Firewall rule '$RULE' deleted"
    else
        echo_warn "Firewall rule '$RULE' does not exist, skipping"
    fi
done

# ============================================================================
# Release static IP (requires region)
# ============================================================================
echo ""
echo_info "Releasing static IP..."

IP_NAME="${RESOURCE_PREFIX}-ip"
if [[ "$ZONE_AVAILABLE" == "true" ]]; then
    if gcloud compute addresses describe "$IP_NAME" --region="$GCP_REGION" &>/dev/null; then
        gcloud compute addresses delete "$IP_NAME" --region="$GCP_REGION" --quiet || echo_warn "Failed to release IP (may be in use)"
        echo_info "Static IP '$IP_NAME' released"
    else
        echo_warn "Static IP '$IP_NAME' does not exist in region $GCP_REGION, skipping"
    fi
else
    # Try to find IP in any region
    local_ip_region=$(gcloud compute addresses list --project="$GCP_PROJECT_ID" \
        --filter="name=${IP_NAME}" --format="value(region)" --quiet 2>/dev/null | head -1 | rev | cut -d'/' -f1 | rev)
    if [[ -n "$local_ip_region" ]]; then
        gcloud compute addresses delete "$IP_NAME" --region="$local_ip_region" --quiet || echo_warn "Failed to release IP"
        echo_info "Static IP '$IP_NAME' released from region $local_ip_region"
    else
        echo_warn "Static IP '$IP_NAME' not found in any region, skipping"
    fi
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=============================================="
echo_info "Infrastructure Teardown Complete!"
echo "=============================================="
echo ""
echo "Removed resources:"
echo "  - Persistent disk: ${RESOURCE_PREFIX}-data"
echo "  - Firewall rules:  ${RESOURCE_PREFIX}-allow-http, ${RESOURCE_PREFIX}-allow-ssh"
echo "  - Static IP:       ${RESOURCE_PREFIX}-ip"
echo ""
