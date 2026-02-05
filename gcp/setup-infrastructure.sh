#!/bin/bash
# Setup GCP infrastructure for qm-nmr-calc Spot VM deployment
#
# Creates:
#   - Static external IP (for DNS configuration)
#   - Firewall rules (HTTP, HTTPS, SSH)
#   - Persistent disk (for job data and Let's Encrypt certs)
#
# Prerequisites:
#   - gcloud CLI installed and authenticated (gcloud auth login)
#   - config.sh created from config.sh.example with your project ID
#
# Usage:
#   cp config.sh.example config.sh
#   # Edit config.sh with your GCP_PROJECT_ID
#   ./setup-infrastructure.sh

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

# Source configuration
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

echo_info "Using configuration:"
echo "  Project:  $GCP_PROJECT_ID"
echo "  Region:   $GCP_REGION"
echo "  Zone:     $GCP_ZONE"
echo "  Prefix:   $RESOURCE_PREFIX"
echo "  Disk:     ${DISK_SIZE_GB}GB SSD"
echo ""

# Check gcloud authentication
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
echo_info "Setting project to $GCP_PROJECT_ID..."
gcloud config set project "$GCP_PROJECT_ID" --quiet

# ============================================================================
# Check Compute Engine API is enabled
# ============================================================================
echo_info "Checking Compute Engine API..."
COMPUTE_API=$(gcloud services list --enabled --filter="name:compute.googleapis.com" --format="value(name)" 2>/dev/null || true)
if [[ -z "$COMPUTE_API" ]]; then
    echo_warn "Compute Engine API is not enabled. Enabling now..."
    if gcloud services enable compute.googleapis.com; then
        echo_info "Compute Engine API enabled successfully"
        sleep 5  # Wait for API to become active
    else
        echo_error "Failed to enable Compute Engine API"
        echo "Enable manually at: https://console.cloud.google.com/apis/library/compute.googleapis.com?project=$GCP_PROJECT_ID"
        exit 1
    fi
else
    echo_info "Compute Engine API is enabled"
fi

# ============================================================================
# INFRA-01: Reserve static external IP
# ============================================================================
echo ""
echo_info "=== INFRA-01: Static External IP ==="

IP_NAME="${RESOURCE_PREFIX}-ip"
EXISTING_IP=$(gcloud compute addresses describe "$IP_NAME" --region="$GCP_REGION" --format="value(address)" 2>/dev/null || true)

if [[ -n "$EXISTING_IP" ]]; then
    echo_warn "Static IP '$IP_NAME' already exists: $EXISTING_IP"
else
    echo_info "Creating static IP '$IP_NAME' in region $GCP_REGION..."
    gcloud compute addresses create "$IP_NAME" --region="$GCP_REGION"
    EXISTING_IP=$(gcloud compute addresses describe "$IP_NAME" --region="$GCP_REGION" --format="value(address)")
    echo_info "Static IP created: $EXISTING_IP"
fi

STATIC_IP="$EXISTING_IP"

# ============================================================================
# INFRA-02: Create firewall rules
# ============================================================================
echo ""
echo_info "=== INFRA-02: Firewall Rules ==="

VM_TAG="${RESOURCE_PREFIX}-vm"

# HTTP (port 80)
HTTP_RULE="${RESOURCE_PREFIX}-allow-http"
if gcloud compute firewall-rules describe "$HTTP_RULE" &>/dev/null; then
    echo_warn "Firewall rule '$HTTP_RULE' already exists"
else
    echo_info "Creating firewall rule '$HTTP_RULE' (TCP:80)..."
    gcloud compute firewall-rules create "$HTTP_RULE" \
        --direction=INGRESS \
        --priority=1000 \
        --network=default \
        --action=ALLOW \
        --rules=tcp:80 \
        --source-ranges=0.0.0.0/0 \
        --target-tags="$VM_TAG"
    echo_info "Firewall rule '$HTTP_RULE' created"
fi

# HTTPS (port 443)
HTTPS_RULE="${RESOURCE_PREFIX}-allow-https"
if gcloud compute firewall-rules describe "$HTTPS_RULE" &>/dev/null; then
    echo_warn "Firewall rule '$HTTPS_RULE' already exists"
else
    echo_info "Creating firewall rule '$HTTPS_RULE' (TCP:443)..."
    gcloud compute firewall-rules create "$HTTPS_RULE" \
        --direction=INGRESS \
        --priority=1000 \
        --network=default \
        --action=ALLOW \
        --rules=tcp:443 \
        --source-ranges=0.0.0.0/0 \
        --target-tags="$VM_TAG"
    echo_info "Firewall rule '$HTTPS_RULE' created"
fi

# SSH (port 22)
SSH_RULE="${RESOURCE_PREFIX}-allow-ssh"
if gcloud compute firewall-rules describe "$SSH_RULE" &>/dev/null; then
    echo_warn "Firewall rule '$SSH_RULE' already exists"
else
    echo_info "Creating firewall rule '$SSH_RULE' (TCP:22)..."
    gcloud compute firewall-rules create "$SSH_RULE" \
        --direction=INGRESS \
        --priority=1000 \
        --network=default \
        --action=ALLOW \
        --rules=tcp:22 \
        --source-ranges=0.0.0.0/0 \
        --target-tags="$VM_TAG"
    echo_info "Firewall rule '$SSH_RULE' created"
fi

# ============================================================================
# INFRA-03: Create persistent disk
# ============================================================================
echo ""
echo_info "=== INFRA-03: Persistent Disk ==="

DISK_NAME="${RESOURCE_PREFIX}-data"
if gcloud compute disks describe "$DISK_NAME" --zone="$GCP_ZONE" &>/dev/null; then
    echo_warn "Persistent disk '$DISK_NAME' already exists"
    DISK_SIZE=$(gcloud compute disks describe "$DISK_NAME" --zone="$GCP_ZONE" --format="value(sizeGb)")
    echo_info "Existing disk size: ${DISK_SIZE}GB"
else
    echo_info "Creating persistent SSD disk '$DISK_NAME' (${DISK_SIZE_GB}GB)..."
    gcloud compute disks create "$DISK_NAME" \
        --zone="$GCP_ZONE" \
        --type=pd-ssd \
        --size="${DISK_SIZE_GB}GB"
    echo_info "Persistent disk '$DISK_NAME' created"
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=============================================="
echo_info "Infrastructure Setup Complete!"
echo "=============================================="
echo ""
echo "Resources created:"
echo "  Static IP:       $STATIC_IP"
echo "  Firewall rules:  $HTTP_RULE, $HTTPS_RULE, $SSH_RULE"
echo "  Persistent disk: $DISK_NAME (${DISK_SIZE_GB}GB SSD)"
echo ""
echo "Next steps:"
echo "  1. Configure DNS: Point your domain A record to $STATIC_IP"
echo "  2. Run Phase 46 (VM Deployment) to create the Spot VM"
echo ""
echo "To verify resources:"
echo "  gcloud compute addresses list --filter=\"name:$IP_NAME\""
echo "  gcloud compute firewall-rules list --filter=\"name~$RESOURCE_PREFIX\""
echo "  gcloud compute disks list --filter=\"name:$DISK_NAME\""
echo ""
