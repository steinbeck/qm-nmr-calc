#!/bin/bash
# Check DNS propagation for qm-nmr-calc deployment
#
# Verifies that your domain resolves to the expected static IP address.
# Run this after configuring DNS and before deploying the VM.
#
# Usage:
#   ./check-dns.sh                    # Uses domain from last deploy
#   ./check-dns.sh nmr.example.com    # Check specific domain

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

# ============================================================================
# Load configuration
# ============================================================================
if [[ ! -f "config.sh" ]]; then
    echo_error "config.sh not found"
    echo "Run ./check-prerequisites.sh first"
    exit 1
fi

source ./config.sh

# ============================================================================
# Get domain from argument or prompt
# ============================================================================
DOMAIN="${1:-}"

if [[ -z "$DOMAIN" ]]; then
    echo ""
    echo "=============================================="
    echo "  DNS Propagation Check"
    echo "=============================================="
    echo ""
    read -r -p "Enter domain to check: " DOMAIN
fi

if [[ -z "$DOMAIN" ]]; then
    echo_error "No domain provided"
    exit 1
fi

# ============================================================================
# Get expected static IP
# ============================================================================
echo ""
echo_check "Looking up static IP..."

# Set project
gcloud config set project "$GCP_PROJECT_ID" --quiet 2>/dev/null || true

IP_NAME="${RESOURCE_PREFIX}-ip"
EXPECTED_IP=$(gcloud compute addresses describe "$IP_NAME" \
    --region="$GCP_REGION" \
    --format="value(address)" 2>/dev/null || true)

if [[ -z "$EXPECTED_IP" ]]; then
    echo_error "Static IP '$IP_NAME' not found"
    echo ""
    echo "Run ./setup-infrastructure.sh first to create the static IP"
    exit 1
fi

echo_info "Expected IP: $EXPECTED_IP"

# ============================================================================
# Check DNS resolution
# ============================================================================
echo ""
echo_check "Checking DNS resolution for $DOMAIN..."

# Try multiple DNS methods for reliability
RESOLVED_IP=""

# Method 1: dig (most reliable)
if command -v dig &>/dev/null; then
    RESOLVED_IP=$(dig +short "$DOMAIN" A 2>/dev/null | head -1 || true)
fi

# Method 2: host (fallback)
if [[ -z "$RESOLVED_IP" ]] && command -v host &>/dev/null; then
    RESOLVED_IP=$(host "$DOMAIN" 2>/dev/null | grep "has address" | head -1 | awk '{print $4}' || true)
fi

# Method 3: nslookup (fallback)
if [[ -z "$RESOLVED_IP" ]] && command -v nslookup &>/dev/null; then
    RESOLVED_IP=$(nslookup "$DOMAIN" 2>/dev/null | grep "Address:" | tail -1 | awk '{print $2}' || true)
fi

# Method 4: getent (Linux)
if [[ -z "$RESOLVED_IP" ]] && command -v getent &>/dev/null; then
    RESOLVED_IP=$(getent hosts "$DOMAIN" 2>/dev/null | awk '{print $1}' || true)
fi

# ============================================================================
# Report results
# ============================================================================
echo ""
echo "=============================================="
echo "  DNS Check Results"
echo "=============================================="
echo ""
echo "Domain:      $DOMAIN"
echo "Expected IP: $EXPECTED_IP"

if [[ -z "$RESOLVED_IP" ]]; then
    echo -e "Resolved IP: ${RED}NOT FOUND${NC}"
    echo ""
    echo_error "Domain does not resolve"
    echo ""
    echo "Possible causes:"
    echo "  - DNS record not created yet"
    echo "  - DNS propagation not complete (wait 5-15 minutes)"
    echo "  - Typo in domain name"
    echo ""
    echo "Required DNS record:"
    echo "  Type: A"
    echo "  Name: $DOMAIN"
    echo "  Value: $EXPECTED_IP"
    echo ""
    echo "Check propagation at: https://dnschecker.org/#A/$DOMAIN"
    exit 1
elif [[ "$RESOLVED_IP" == "$EXPECTED_IP" ]]; then
    echo -e "Resolved IP: ${GREEN}$RESOLVED_IP${NC}"
    echo ""
    echo_info "DNS is correctly configured!"
    echo ""
    echo "You can now deploy the VM:"
    echo "  ./deploy-vm.sh"
    echo ""
else
    echo -e "Resolved IP: ${YELLOW}$RESOLVED_IP${NC}"
    echo ""
    echo_warn "Domain resolves to a different IP address"
    echo ""
    echo "This could mean:"
    echo "  - DNS record points to an old IP"
    echo "  - Cloudflare proxy is enabled (should be disabled for initial setup)"
    echo "  - DNS propagation still in progress"
    echo ""
    echo "Update your DNS A record to point to: $EXPECTED_IP"
    echo ""
    echo "Check propagation at: https://dnschecker.org/#A/$DOMAIN"
    exit 1
fi
