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

# Source configuration
if [[ ! -f "config.sh" ]]; then
    echo_error "config.sh not found!"
    echo "Copy config.sh.example to config.sh and set your GCP_PROJECT_ID:"
    echo "  cp config.sh.example config.sh"
    exit 1
fi

source ./config.sh

# Validate required variables
if [[ -z "${GCP_PROJECT_ID:-}" || "$GCP_PROJECT_ID" == "your-project-id" ]]; then
    echo_error "GCP_PROJECT_ID is not set or still has default value"
    exit 1
fi

echo_info "Teardown configuration:"
echo "  Project:  $GCP_PROJECT_ID"
echo "  Region:   $GCP_REGION"
echo "  Zone:     $GCP_ZONE"
echo "  Prefix:   $RESOURCE_PREFIX"
echo ""

# Set project
gcloud config set project "$GCP_PROJECT_ID" --quiet

# Confirmation prompt
if [[ "$FORCE" != "true" ]]; then
    echo ""
    echo_warn "WARNING: This will delete the following resources:"
    echo "  - Persistent disk '${RESOURCE_PREFIX}-data' with ALL job data"
    echo "  - Firewall rules for HTTP, HTTPS, SSH"
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
# Delete persistent disk (first, as it's the most important)
# ============================================================================
DISK_NAME="${RESOURCE_PREFIX}-data"
echo ""
echo_info "Deleting persistent disk '$DISK_NAME'..."
if gcloud compute disks describe "$DISK_NAME" --zone="$GCP_ZONE" &>/dev/null; then
    gcloud compute disks delete "$DISK_NAME" --zone="$GCP_ZONE" --quiet || echo_warn "Failed to delete disk (may be attached to VM)"
    echo_info "Persistent disk '$DISK_NAME' deleted"
else
    echo_warn "Persistent disk '$DISK_NAME' does not exist, skipping"
fi

# ============================================================================
# Delete firewall rules
# ============================================================================
echo ""
echo_info "Deleting firewall rules..."

HTTP_RULE="${RESOURCE_PREFIX}-allow-http"
HTTPS_RULE="${RESOURCE_PREFIX}-allow-https"
SSH_RULE="${RESOURCE_PREFIX}-allow-ssh"

for RULE in "$HTTP_RULE" "$HTTPS_RULE" "$SSH_RULE"; do
    if gcloud compute firewall-rules describe "$RULE" &>/dev/null; then
        gcloud compute firewall-rules delete "$RULE" --quiet || echo_warn "Failed to delete firewall rule '$RULE'"
        echo_info "Firewall rule '$RULE' deleted"
    else
        echo_warn "Firewall rule '$RULE' does not exist, skipping"
    fi
done

# ============================================================================
# Release static IP
# ============================================================================
echo ""
echo_info "Releasing static IP..."

IP_NAME="${RESOURCE_PREFIX}-ip"
if gcloud compute addresses describe "$IP_NAME" --region="$GCP_REGION" &>/dev/null; then
    gcloud compute addresses delete "$IP_NAME" --region="$GCP_REGION" --quiet || echo_warn "Failed to release IP (may be in use)"
    echo_info "Static IP '$IP_NAME' released"
else
    echo_warn "Static IP '$IP_NAME' does not exist, skipping"
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
echo "  - Persistent disk: $DISK_NAME"
echo "  - Firewall rules:  $HTTP_RULE, $HTTPS_RULE, $SSH_RULE"
echo "  - Static IP:       $IP_NAME"
echo ""
echo "Note: Remember to update your DNS records if you had configured them"
echo ""
