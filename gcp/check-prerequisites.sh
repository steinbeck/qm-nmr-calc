#!/bin/bash
# Check GCP prerequisites for qm-nmr-calc deployment
#
# Verifies all prerequisites are met before attempting deployment:
#   - gcloud CLI installed and authenticated
#   - Project exists and is accessible
#   - Billing enabled on project
#   - Compute Engine API enabled
#   - Required permissions
#
# Usage:
#   ./check-prerequisites.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo_info() { echo -e "${GREEN}[✓]${NC} $1"; }
echo_warn() { echo -e "${YELLOW}[!]${NC} $1"; }
echo_error() { echo -e "${RED}[✗]${NC} $1"; }
echo_check() { echo -e "${BLUE}[?]${NC} $1"; }

ERRORS=0
WARNINGS=0

echo ""
echo "=============================================="
echo "  GCP Prerequisites Check"
echo "=============================================="
echo ""

# ============================================================================
# Check 1: gcloud CLI installed
# ============================================================================
echo_check "Checking gcloud CLI installation..."
if ! command -v gcloud &>/dev/null; then
    echo_error "gcloud CLI is not installed"
    echo ""
    echo "  Install the Google Cloud SDK:"
    echo "    macOS:   brew install google-cloud-sdk"
    echo "    Linux:   curl https://sdk.cloud.google.com | bash"
    echo "    Windows: https://cloud.google.com/sdk/docs/install"
    echo ""
    ((ERRORS++))
else
    GCLOUD_VERSION=$(gcloud --version 2>/dev/null | head -1)
    echo_info "gcloud CLI installed: $GCLOUD_VERSION"
fi

# ============================================================================
# Check 2: gcloud authenticated
# ============================================================================
echo_check "Checking gcloud authentication..."
ACTIVE_ACCOUNT=$(gcloud auth list --filter=status:ACTIVE --format="value(account)" 2>/dev/null || true)
if [[ -z "$ACTIVE_ACCOUNT" ]]; then
    echo_error "Not authenticated with gcloud"
    echo ""
    echo "  Run: gcloud auth login"
    echo ""
    ((ERRORS++))
else
    echo_info "Authenticated as: $ACTIVE_ACCOUNT"
fi

# ============================================================================
# Check 3: config.sh exists
# ============================================================================
echo_check "Checking config.sh..."
if [[ ! -f "config.sh" ]]; then
    echo_error "config.sh not found"
    echo ""
    echo "  Create config.sh from template:"
    echo "    cp config.sh.example config.sh"
    echo "    # Edit config.sh and set GCP_PROJECT_ID"
    echo ""
    ((ERRORS++))
else
    source ./config.sh

    # Check GCP_PROJECT_ID is set
    if [[ -z "${GCP_PROJECT_ID:-}" || "$GCP_PROJECT_ID" == "your-project-id" ]]; then
        echo_error "GCP_PROJECT_ID not configured in config.sh"
        echo ""
        echo "  Edit config.sh and set your actual GCP project ID"
        echo ""
        ((ERRORS++))
    else
        echo_info "config.sh loaded: project=$GCP_PROJECT_ID"
    fi
fi

# Exit early if critical errors
if [[ $ERRORS -gt 0 && (-z "${GCP_PROJECT_ID:-}" || "$GCP_PROJECT_ID" == "your-project-id") ]]; then
    echo ""
    echo "=============================================="
    echo_error "Prerequisites check failed"
    echo "=============================================="
    echo ""
    echo "Fix the errors above and run this script again."
    echo ""
    exit 1
fi

# Set project for remaining checks
gcloud config set project "$GCP_PROJECT_ID" --quiet 2>/dev/null || true

# ============================================================================
# Check 4: Project exists and is accessible
# ============================================================================
echo_check "Checking project access..."
if ! gcloud projects describe "$GCP_PROJECT_ID" &>/dev/null; then
    echo_error "Cannot access project '$GCP_PROJECT_ID'"
    echo ""
    echo "  Possible causes:"
    echo "    - Project does not exist"
    echo "    - You don't have access to the project"
    echo "    - Project ID is incorrect"
    echo ""
    echo "  To list your accessible projects:"
    echo "    gcloud projects list"
    echo ""
    ((ERRORS++))
else
    PROJECT_NAME=$(gcloud projects describe "$GCP_PROJECT_ID" --format="value(name)" 2>/dev/null || echo "unknown")
    echo_info "Project accessible: $PROJECT_NAME ($GCP_PROJECT_ID)"
fi

# ============================================================================
# Check 5: Billing enabled
# ============================================================================
echo_check "Checking billing status..."
BILLING_ACCOUNT=$(gcloud billing projects describe "$GCP_PROJECT_ID" --format="value(billingAccountName)" 2>/dev/null || true)
if [[ -z "$BILLING_ACCOUNT" ]]; then
    echo_error "Billing is not enabled for project '$GCP_PROJECT_ID'"
    echo ""
    echo "  Enable billing at:"
    echo "    https://console.cloud.google.com/billing/linkedaccount?project=$GCP_PROJECT_ID"
    echo ""
    ((ERRORS++))
else
    echo_info "Billing enabled"
fi

# ============================================================================
# Check 6: Compute Engine API enabled
# ============================================================================
echo_check "Checking Compute Engine API..."
COMPUTE_API=$(gcloud services list --enabled --filter="name:compute.googleapis.com" --format="value(name)" 2>/dev/null || true)
if [[ -z "$COMPUTE_API" ]]; then
    echo_error "Compute Engine API is not enabled"
    echo ""
    echo "  Enable it with:"
    echo "    gcloud services enable compute.googleapis.com"
    echo ""
    echo "  Or via console:"
    echo "    https://console.cloud.google.com/apis/library/compute.googleapis.com?project=$GCP_PROJECT_ID"
    echo ""
    ((ERRORS++))
else
    echo_info "Compute Engine API enabled"
fi

# ============================================================================
# Check 7: Region/zone validity (optional check)
# ============================================================================
if [[ -n "${GCP_REGION:-}" && -n "${GCP_ZONE:-}" ]]; then
    echo_check "Checking region/zone validity..."
    if gcloud compute zones describe "$GCP_ZONE" &>/dev/null; then
        echo_info "Zone valid: $GCP_ZONE"
    else
        echo_warn "Zone '$GCP_ZONE' may not exist or be available"
        echo "  List available zones: gcloud compute zones list"
        ((WARNINGS++))
    fi
fi

# ============================================================================
# Check 8: Sufficient quota (basic check)
# ============================================================================
echo_check "Checking compute quotas..."
# Just check if we can query quotas - detailed check would require parsing
if gcloud compute regions describe "${GCP_REGION:-us-central1}" --format="value(quotas)" &>/dev/null; then
    echo_info "Quota information accessible (run ./deploy-vm.sh to verify specific requirements)"
else
    echo_warn "Could not verify quotas (may need quota increases for larger VMs)"
    ((WARNINGS++))
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=============================================="
if [[ $ERRORS -eq 0 ]]; then
    echo_info "All prerequisites passed!"
    echo "=============================================="
    echo ""
    echo "You're ready to deploy. Run these commands in order:"
    echo ""
    echo "  1. ./setup-infrastructure.sh   # Create static IP, firewall, disk"
    echo "  2. Configure DNS (point domain to static IP shown above)"
    echo "  3. ./deploy-vm.sh              # Create and start the VM"
    echo ""
    if [[ $WARNINGS -gt 0 ]]; then
        echo "Note: $WARNINGS warning(s) - review them above if deployment fails."
        echo ""
    fi
else
    echo_error "Prerequisites check failed: $ERRORS error(s)"
    echo "=============================================="
    echo ""
    echo "Fix the errors above and run this script again:"
    echo "  ./check-prerequisites.sh"
    echo ""
    exit 1
fi
