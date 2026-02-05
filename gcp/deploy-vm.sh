#!/bin/bash
# Deploy qm-nmr-calc to GCP Spot VM
#
# This script creates a cost-optimized Spot VM with:
#   - Interactive prompts for region, zone, machine type, and domain
#   - Cost estimation before creation
#   - DNS validation before deployment
#   - Automatic Docker installation and container deployment via startup script
#   - Graceful shutdown handling during preemption
#
# Prerequisites:
#   - gcloud CLI installed and authenticated (gcloud auth login)
#   - config.sh created from config.sh.example with your project ID
#   - Infrastructure created (run ./setup-infrastructure.sh first)
#   - DNS configured (domain pointing to static IP)
#
# Usage:
#   ./deploy-vm.sh
#   ./deploy-vm.sh --skip-dns-check    # Skip DNS validation (not recommended)
#
# First time? Run these in order:
#   1. ./check-prerequisites.sh
#   2. ./setup-infrastructure.sh
#   3. Configure DNS (point domain to static IP)
#   4. ./check-dns.sh your-domain.com
#   5. ./deploy-vm.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Parse arguments
SKIP_DNS_CHECK=false
for arg in "$@"; do
    case $arg in
        --skip-dns-check)
            SKIP_DNS_CHECK=true
            shift
            ;;
    esac
done

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
echo_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
echo_error() { echo -e "${RED}[ERROR]${NC} $1"; }
echo_prompt() { echo -e "${BLUE}[INPUT]${NC} $1"; }

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
# Check Compute Engine API is enabled
# ============================================================================
echo_info "Checking Compute Engine API..."
COMPUTE_API=$(gcloud services list --enabled --filter="name:compute.googleapis.com" --format="value(name)" 2>/dev/null || true)
if [[ -z "$COMPUTE_API" ]]; then
    echo_error "Compute Engine API is not enabled"
    echo ""
    echo "Enable it with:"
    echo "  gcloud services enable compute.googleapis.com"
    echo ""
    echo "Or run: ./check-prerequisites.sh"
    exit 1
fi

# ============================================================================
# Interactive prompts
# ============================================================================
echo ""
echo "=============================================="
echo "  GCP Spot VM Deployment for qm-nmr-calc"
echo "=============================================="
echo ""

# Region selection
echo_prompt "Select region"
read -r -p "Region [${GCP_REGION}]: " INPUT_REGION
SELECTED_REGION="${INPUT_REGION:-${GCP_REGION}}"

# Zone selection
echo_prompt "Select zone"
read -r -p "Zone [${GCP_ZONE}]: " INPUT_ZONE
SELECTED_ZONE="${INPUT_ZONE:-${GCP_ZONE}}"

# Machine type selection with cost information
echo ""
echo "Available machine types (Spot pricing ~60-91% discount):"
echo "  e2-standard-2   (2 vCPU, 8 GB)   ~\$15-25/month  - Light workloads"
echo "  e2-standard-4   (4 vCPU, 16 GB)  ~\$30-50/month  - Recommended"
echo "  n2-standard-4   (4 vCPU, 16 GB)  ~\$50-80/month  - Higher performance"
echo "  c2-standard-4   (4 vCPU, 16 GB)  ~\$70-100/month - Compute-optimized"
echo ""
echo_prompt "Select machine type"
read -r -p "Machine type [e2-standard-4]: " INPUT_MACHINE
MACHINE_TYPE="${INPUT_MACHINE:-e2-standard-4}"

# Domain prompt (required for HTTPS)
echo ""
echo_prompt "Enter domain for HTTPS (e.g., nmr.example.com)"
echo "This domain must point to your static IP for Let's Encrypt to work."
read -r -p "Domain: " DOMAIN

# Validate DOMAIN is provided
if [[ -z "$DOMAIN" ]]; then
    echo_error "Domain is required for HTTPS deployment"
    exit 1
fi

# ============================================================================
# Validate DNS configuration
# ============================================================================
if [[ "$SKIP_DNS_CHECK" != "true" ]]; then
    echo ""
    echo_info "Validating DNS configuration..."

    # Get expected static IP
    IP_NAME="${RESOURCE_PREFIX}-ip"
    EXPECTED_IP=$(gcloud compute addresses describe "$IP_NAME" \
        --region="$SELECTED_REGION" \
        --format="value(address)" 2>/dev/null || true)

    if [[ -z "$EXPECTED_IP" ]]; then
        echo_error "Static IP not found. Run ./setup-infrastructure.sh first."
        exit 1
    fi

    # Check DNS resolution
    RESOLVED_IP=""
    if command -v dig &>/dev/null; then
        RESOLVED_IP=$(dig +short "$DOMAIN" A 2>/dev/null | head -1 || true)
    elif command -v host &>/dev/null; then
        RESOLVED_IP=$(host "$DOMAIN" 2>/dev/null | grep "has address" | head -1 | awk '{print $4}' || true)
    elif command -v nslookup &>/dev/null; then
        RESOLVED_IP=$(nslookup "$DOMAIN" 2>/dev/null | grep "Address:" | tail -1 | awk '{print $2}' || true)
    fi

    if [[ -z "$RESOLVED_IP" ]]; then
        echo_error "Domain '$DOMAIN' does not resolve"
        echo ""
        echo "Configure your DNS A record to point to: $EXPECTED_IP"
        echo "Then wait for propagation (5-15 minutes) and try again."
        echo ""
        echo "Check propagation at: https://dnschecker.org/#A/$DOMAIN"
        echo ""
        echo "To skip this check (not recommended):"
        echo "  ./deploy-vm.sh --skip-dns-check"
        exit 1
    elif [[ "$RESOLVED_IP" != "$EXPECTED_IP" ]]; then
        echo_warn "Domain '$DOMAIN' resolves to $RESOLVED_IP (expected: $EXPECTED_IP)"
        echo ""
        echo "This may cause Let's Encrypt certificate issues."
        echo ""
        read -r -p "Continue anyway? (yes/no) [no]: " DNS_CONTINUE
        DNS_CONTINUE="${DNS_CONTINUE:-no}"
        if [[ "$DNS_CONTINUE" != "yes" ]]; then
            echo_info "Deployment cancelled. Fix DNS and try again."
            exit 0
        fi
    else
        echo_info "DNS verified: $DOMAIN -> $EXPECTED_IP"
    fi
fi

# ============================================================================
# Cost estimation and confirmation
# ============================================================================
echo ""
echo "=============================================="
echo "  Configuration Summary"
echo "=============================================="
echo ""
echo "Project:      $GCP_PROJECT_ID"
echo "Region:       $SELECTED_REGION"
echo "Zone:         $SELECTED_ZONE"
echo "Machine type: $MACHINE_TYPE"
echo "Domain:       $DOMAIN"
echo ""
echo "Estimated monthly cost (Spot VM, ~730 hours uptime):"
case "$MACHINE_TYPE" in
    e2-standard-2)
        echo "  ~\$15-25/month"
        ;;
    e2-standard-4)
        echo "  ~\$30-50/month"
        ;;
    n2-standard-4)
        echo "  ~\$50-80/month"
        ;;
    c2-standard-4)
        echo "  ~\$70-100/month"
        ;;
    *)
        echo "  Cost varies by machine type"
        ;;
esac
echo ""
echo "Note: Spot pricing is dynamic and can vary. Actual costs depend on"
echo "market conditions and uptime. VMs may be preempted with 30s notice."
echo "Visit https://cloud.google.com/compute/vm-instance-pricing for details."
echo ""

# Confirmation prompt
read -r -p "Continue with VM creation? (yes/no) [yes]: " CONFIRM
CONFIRM="${CONFIRM:-yes}"

if [[ "$CONFIRM" != "yes" ]]; then
    echo_info "Deployment cancelled by user"
    exit 0
fi

# ============================================================================
# Get static IP address
# ============================================================================
echo ""
echo_info "Retrieving static IP address..."
IP_NAME="${RESOURCE_PREFIX}-ip"
STATIC_IP=$(gcloud compute addresses describe "$IP_NAME" \
    --region="$SELECTED_REGION" \
    --format="value(address)" 2>/dev/null || true)

if [[ -z "$STATIC_IP" ]]; then
    echo_error "Static IP '$IP_NAME' not found in region $SELECTED_REGION"
    echo "Run ./setup-infrastructure.sh first to create required infrastructure"
    exit 1
fi

echo_info "Static IP: $STATIC_IP"

# ============================================================================
# Check if VM already exists
# ============================================================================
VM_NAME="${RESOURCE_PREFIX}-vm"
if gcloud compute instances describe "$VM_NAME" --zone="$SELECTED_ZONE" &>/dev/null; then
    echo_warn "VM '$VM_NAME' already exists in zone $SELECTED_ZONE"
    echo ""
    read -r -p "Delete and recreate VM? (yes/no) [no]: " RECREATE
    RECREATE="${RECREATE:-no}"

    if [[ "$RECREATE" == "yes" ]]; then
        echo_info "Deleting existing VM..."
        gcloud compute instances delete "$VM_NAME" --zone="$SELECTED_ZONE" --quiet
        echo_info "VM deleted"
    else
        echo_info "Deployment cancelled - VM already exists"
        exit 0
    fi
fi

# ============================================================================
# Verify required scripts exist
# ============================================================================
if [[ ! -f "startup.sh" || ! -f "shutdown.sh" ]]; then
    echo_error "Required scripts not found (startup.sh, shutdown.sh)"
    echo "These should be in the gcp/ directory"
    exit 1
fi

# ============================================================================
# Create Spot VM
# ============================================================================
echo ""
echo_info "Creating Spot VM '$VM_NAME'..."
echo_info "This will take 2-3 minutes..."

DISK_NAME="${RESOURCE_PREFIX}-data"

gcloud compute instances create "$VM_NAME" \
    --zone="$SELECTED_ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --provisioning-model=SPOT \
    --instance-termination-action=STOP \
    --tags="${RESOURCE_PREFIX}-vm" \
    --network-interface="address=$STATIC_IP" \
    --boot-disk-size=20GB \
    --boot-disk-type=pd-balanced \
    --image-family=debian-12 \
    --image-project=debian-cloud \
    --disk="name=${DISK_NAME},mode=rw,device-name=data-disk,boot=no,auto-delete=no" \
    --metadata-from-file=startup-script=startup.sh,shutdown-script=shutdown.sh \
    --metadata="DOMAIN=${DOMAIN}"

echo_info "VM created successfully!"

# ============================================================================
# Post-creation summary
# ============================================================================
echo ""
echo "=============================================="
echo "  VM Deployment Complete!"
echo "=============================================="
echo ""
echo "VM Name:    $VM_NAME"
echo "Zone:       $SELECTED_ZONE"
echo "Static IP:  $STATIC_IP"
echo "Domain:     $DOMAIN"
echo ""
echo "Next steps:"
echo ""
echo "1. Configure DNS:"
echo "   Point your domain A record to $STATIC_IP"
echo "   Example DNS record:"
echo "     $DOMAIN.  A  $STATIC_IP"
echo ""
echo "2. Wait for startup (2-3 minutes):"
echo "   The VM is installing Docker and starting containers."
echo "   Monitor startup progress:"
echo "     gcloud compute ssh $VM_NAME --zone=$SELECTED_ZONE --command=\"tail -f /var/log/startup-script.log\""
echo ""
echo "3. Verify deployment:"
echo "   Check running containers:"
echo "     gcloud compute ssh $VM_NAME --zone=$SELECTED_ZONE --command=\"docker compose -f /opt/qm-nmr-calc/docker-compose.yml ps\""
echo ""
echo "4. Access your deployment:"
echo "   https://$DOMAIN"
echo "   (Let's Encrypt certificate will be issued automatically once DNS propagates)"
echo ""
echo "To SSH into the VM:"
echo "  gcloud compute ssh $VM_NAME --zone=$SELECTED_ZONE"
echo ""
echo "To view logs:"
echo "  gcloud compute ssh $VM_NAME --zone=$SELECTED_ZONE --command=\"docker compose -f /opt/qm-nmr-calc/docker-compose.yml logs -f\""
echo ""
