#!/bin/bash
# Verify GCP deployment for qm-nmr-calc
#
# Checks that all components are running correctly:
#   - VM is running
#   - Containers are healthy
#   - HTTPS is working
#   - API responds correctly
#
# Usage:
#   ./verify-deployment.sh

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

# ============================================================================
# Load configuration
# ============================================================================
if [[ ! -f "config.sh" ]]; then
    echo_error "config.sh not found"
    exit 1
fi

source ./config.sh

echo ""
echo "=============================================="
echo "  Deployment Verification"
echo "=============================================="
echo ""

# Set project
gcloud config set project "$GCP_PROJECT_ID" --quiet 2>/dev/null || true

# ============================================================================
# Check 1: VM exists and is running
# ============================================================================
VM_NAME="${RESOURCE_PREFIX}-vm"
echo_check "Checking VM status..."

if ! gcloud compute instances describe "$VM_NAME" --zone="$GCP_ZONE" &>/dev/null; then
    echo_error "VM '$VM_NAME' not found"
    echo "  Run ./deploy-vm.sh to create the VM"
    exit 1
fi

VM_STATUS=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(status)")

if [[ "$VM_STATUS" == "RUNNING" ]]; then
    echo_info "VM is running"
else
    echo_error "VM is not running (status: $VM_STATUS)"
    echo "  Run ./start-vm.sh to start the VM"
    ((ERRORS++))
fi

# Get VM details
EXTERNAL_IP=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(networkInterfaces[0].accessConfigs[0].natIP)" 2>/dev/null || echo "")

DOMAIN=$(gcloud compute instances describe "$VM_NAME" \
    --zone="$GCP_ZONE" \
    --format="value(metadata.items[DOMAIN])" 2>/dev/null || echo "")

echo "  External IP: ${EXTERNAL_IP:-N/A}"
echo "  Domain: ${DOMAIN:-N/A}"
echo ""

# ============================================================================
# Check 2: Container status (via SSH)
# ============================================================================
echo_check "Checking container status..."

if [[ "$VM_STATUS" == "RUNNING" ]]; then
    CONTAINER_STATUS=$(gcloud compute ssh "$VM_NAME" --zone="$GCP_ZONE" --quiet \
        --command="docker compose -f /opt/qm-nmr-calc/docker-compose.yml -f /opt/qm-nmr-calc/docker-compose.gcp.yml ps --format '{{.Service}}:{{.Status}}:{{.Health}}'" 2>/dev/null || echo "")

    if [[ -n "$CONTAINER_STATUS" ]]; then
        echo "  Containers:"
        while IFS= read -r line; do
            SERVICE=$(echo "$line" | cut -d: -f1)
            STATUS=$(echo "$line" | cut -d: -f2)
            HEALTH=$(echo "$line" | cut -d: -f3)

            if [[ "$HEALTH" == "healthy" ]] || [[ "$STATUS" == *"Up"* && -z "$HEALTH" ]]; then
                echo -e "    ${GREEN}✓${NC} $SERVICE: $STATUS"
            elif [[ "$STATUS" == *"Up"* ]]; then
                echo -e "    ${YELLOW}!${NC} $SERVICE: $STATUS (health: $HEALTH)"
                ((WARNINGS++))
            else
                echo -e "    ${RED}✗${NC} $SERVICE: $STATUS"
                ((ERRORS++))
            fi
        done <<< "$CONTAINER_STATUS"
    else
        echo_warn "Could not retrieve container status"
        ((WARNINGS++))
    fi
else
    echo_warn "Skipping container check (VM not running)"
fi
echo ""

# ============================================================================
# Check 3: HTTPS endpoint
# ============================================================================
if [[ -n "$DOMAIN" && "$DOMAIN" != "localhost" && "$VM_STATUS" == "RUNNING" ]]; then
    echo_check "Checking HTTPS endpoint..."

    # Try HTTPS first
    HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" --connect-timeout 10 "https://$DOMAIN/health" 2>/dev/null || echo "000")

    if [[ "$HTTP_CODE" == "200" ]]; then
        echo_info "HTTPS endpoint responding (https://$DOMAIN)"
    elif [[ "$HTTP_CODE" == "000" ]]; then
        echo_warn "Could not connect to https://$DOMAIN"
        echo "  This may be normal if:"
        echo "    - DNS hasn't propagated yet (wait 5-15 minutes)"
        echo "    - Let's Encrypt certificate is still being issued"
        echo "    - Firewall rules are not configured"
        ((WARNINGS++))
    else
        echo_warn "HTTPS returned status $HTTP_CODE"
        ((WARNINGS++))
    fi
    echo ""
fi

# ============================================================================
# Check 4: API health (via SSH)
# ============================================================================
if [[ "$VM_STATUS" == "RUNNING" ]]; then
    echo_check "Checking API health..."

    API_HEALTH=$(gcloud compute ssh "$VM_NAME" --zone="$GCP_ZONE" --quiet \
        --command="curl -s http://localhost:8000/health 2>/dev/null || echo 'FAILED'" 2>/dev/null || echo "SSH_FAILED")

    if [[ "$API_HEALTH" == *"healthy"* ]] || [[ "$API_HEALTH" == *"ok"* ]]; then
        echo_info "API health check passed"
    elif [[ "$API_HEALTH" == "SSH_FAILED" ]]; then
        echo_warn "Could not SSH to VM for health check"
        ((WARNINGS++))
    else
        echo_warn "API health check response: $API_HEALTH"
        ((WARNINGS++))
    fi
    echo ""
fi

# ============================================================================
# Check 5: Disk mount
# ============================================================================
if [[ "$VM_STATUS" == "RUNNING" ]]; then
    echo_check "Checking persistent disk mount..."

    DISK_MOUNTED=$(gcloud compute ssh "$VM_NAME" --zone="$GCP_ZONE" --quiet \
        --command="mountpoint -q /mnt/disks/data && echo 'mounted' || echo 'not mounted'" 2>/dev/null || echo "unknown")

    if [[ "$DISK_MOUNTED" == "mounted" ]]; then
        echo_info "Persistent disk mounted at /mnt/disks/data"
    else
        echo_error "Persistent disk not mounted"
        ((ERRORS++))
    fi
    echo ""
fi

# ============================================================================
# Check 6: Startup script completion
# ============================================================================
if [[ "$VM_STATUS" == "RUNNING" ]]; then
    echo_check "Checking startup script completion..."

    STARTUP_LOG=$(gcloud compute ssh "$VM_NAME" --zone="$GCP_ZONE" --quiet \
        --command="tail -5 /var/log/startup-script.log 2>/dev/null || echo 'NO_LOG'" 2>/dev/null || echo "SSH_FAILED")

    if [[ "$STARTUP_LOG" == *"Startup script finished"* ]] || [[ "$STARTUP_LOG" == *"VM Setup Complete"* ]]; then
        echo_info "Startup script completed successfully"
    elif [[ "$STARTUP_LOG" == "NO_LOG" ]] || [[ "$STARTUP_LOG" == "SSH_FAILED" ]]; then
        echo_warn "Could not check startup script log"
        ((WARNINGS++))
    else
        echo_warn "Startup script may still be running"
        echo "  View logs: ./logs-vm.sh"
        ((WARNINGS++))
    fi
    echo ""
fi

# ============================================================================
# Summary
# ============================================================================
echo "=============================================="
if [[ $ERRORS -eq 0 && $WARNINGS -eq 0 ]]; then
    echo_info "All checks passed!"
    echo "=============================================="
    echo ""
    if [[ -n "$DOMAIN" && "$DOMAIN" != "localhost" ]]; then
        echo "Your deployment is ready at:"
        echo "  https://$DOMAIN"
    fi
elif [[ $ERRORS -eq 0 ]]; then
    echo_warn "Deployment OK with $WARNINGS warning(s)"
    echo "=============================================="
    echo ""
    echo "Review warnings above. The deployment may still work."
    if [[ -n "$DOMAIN" && "$DOMAIN" != "localhost" ]]; then
        echo ""
        echo "Try accessing: https://$DOMAIN"
    fi
else
    echo_error "Deployment has $ERRORS error(s) and $WARNINGS warning(s)"
    echo "=============================================="
    echo ""
    echo "Fix the errors above. Useful commands:"
    echo "  ./status-vm.sh   - Check VM status"
    echo "  ./logs-vm.sh     - View container logs"
    echo "  ./ssh-vm.sh      - SSH into VM for debugging"
fi
echo ""
