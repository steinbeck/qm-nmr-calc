#!/bin/bash
# Automated GCP deployment for qm-nmr-calc (v2.7)
#
# Deploys end-to-end from TOML config with zero interactive prompts.
# Replaces v2.6 interactive deploy-vm.sh.
#
# Usage:
#   ./deploy-auto.sh                    # Full deployment
#   ./deploy-auto.sh --dry-run          # Show planned actions without executing
#   ./deploy-auto.sh --config path.toml # Use custom config file
#
# Prerequisites:
#   - gcloud CLI installed and authenticated
#   - config.toml created from config.toml.example
#   - jq installed (apt-get install jq)

set -euo pipefail

# ============================================================================
# Argument parsing
# ============================================================================

DRY_RUN=false
CONFIG_PATH="gcp/config.toml"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            echo "Usage: $0 [--dry-run] [--config <path>]" >&2
            exit 1
            ;;
    esac
done

# Export DRY_RUN so infra.sh sees it
export DRY_RUN

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source all libraries
source "$SCRIPT_DIR/lib/config.sh"
source "$SCRIPT_DIR/lib/pricing.sh"
source "$SCRIPT_DIR/lib/machine.sh"
source "$SCRIPT_DIR/lib/infra.sh"

# Register cleanup trap (from infra.sh)
trap cleanup_on_failure EXIT

# ============================================================================
# Main deployment pipeline (6 steps)
# ============================================================================

main() {
    log_info "=========================================="
    log_info "  qm-nmr-calc Automated GCP Deployment"
    log_info "=========================================="
    [[ "$DRY_RUN" == "true" ]] && log_warn "DRY-RUN MODE: No resources will be created"
    echo ""

    # ── Step 1: Load and validate config ──
    log_info "Step 1/6: Loading configuration..."
    load_config "$CONFIG_PATH" || { log_error "Config validation failed"; exit 1; }
    log_info "  Project:     $GCP_PROJECT_ID"
    log_info "  Prefix:      $RESOURCE_PREFIX"
    log_info "  CPU cores:   $CPU_CORES"
    log_info "  RAM:         ${RAM_GB}GB"
    log_info "  Disk:        ${DISK_SIZE_GB}GB"
    echo ""

    # ── Step 2: Check gcloud authentication ──
    log_info "Step 2/6: Checking gcloud authentication..."
    local active_account
    active_account=$(gcloud auth list --filter=status:ACTIVE --format="value(account)" 2>/dev/null || true)
    if [[ -z "$active_account" ]]; then
        log_error "Not authenticated with gcloud. Run: gcloud auth login"
        exit 1
    fi
    log_info "  Authenticated as: $active_account"
    gcloud config set project "$GCP_PROJECT_ID" --quiet
    echo ""

    # ── Step 3: Select machine and zone ──
    log_info "Step 3/6: Finding cheapest machine type and zone..."
    local machine_json
    machine_json=$(select_machine "$CPU_CORES" "$RAM_GB") || { log_error "Machine selection failed"; exit 1; }

    local machine_type zone region
    machine_type=$(echo "$machine_json" | jq -r '.machine_type')
    zone=$(echo "$machine_json" | jq -r '.zone')
    region=$(echo "$machine_json" | jq -r '.region')

    log_info "  Machine type: $machine_type"
    log_info "  Zone:         $zone"
    log_info "  Region:       $region"
    echo ""

    # Get spot pricing for cost estimate
    local pricing_json spot_price_hourly
    pricing_json=$(get_pricing_table "$CPU_CORES" "$RAM_GB") || true
    spot_price_hourly=$(echo "$pricing_json" | jq -r '.[0].spot_price_hourly // "0.05"')

    # ── Step 4: Display cost estimate ──
    log_info "Step 4/6: Cost estimate..."
    display_cost_estimate "$machine_type" "$spot_price_hourly" "$DISK_SIZE_GB"
    echo ""

    # ── Step 5: Create infrastructure ──
    log_info "Step 5/6: Creating infrastructure..."
    local ip_name="${RESOURCE_PREFIX}-ip"
    local disk_name="${RESOURCE_PREFIX}-data"
    local vm_name="${RESOURCE_PREFIX}-vm"

    local static_ip
    static_ip=$(create_static_ip "$ip_name" "$region")
    create_firewall_rules "$RESOURCE_PREFIX"
    create_persistent_disk "$disk_name" "$zone" "$DISK_SIZE_GB"
    echo ""

    # ── Step 6: Generate startup script and create VM ──
    log_info "Step 6/6: Creating Spot VM..."

    # Generate dynamic startup script using Phase 50 library
    local startup_tmp
    startup_tmp=$(mktemp /tmp/startup-XXXXXX.sh)
    generate_startup "$CPU_CORES" "$RAM_GB" "$RESOURCE_PREFIX" "$DISK_SIZE_GB" > "$startup_tmp"
    log_info "  Startup script generated ($startup_tmp)"

    # Use existing shutdown.sh from gcp/
    local shutdown_script="$SCRIPT_DIR/shutdown.sh"
    if [[ ! -f "$shutdown_script" ]]; then
        log_error "Shutdown script not found: $shutdown_script"
        exit 1
    fi

    create_vm "$vm_name" "$zone" "$machine_type" "$static_ip" "$disk_name" \
        "$startup_tmp" "$shutdown_script" "$RESOURCE_PREFIX"

    # Clean up temp file
    rm -f "$startup_tmp"
    echo ""

    # ── Summary ──
    log_info "=========================================="
    log_info "  Deployment Complete!"
    log_info "=========================================="
    echo ""
    log_info "VM Name:     $vm_name"
    log_info "Zone:        $zone"
    log_info "Machine:     $machine_type"
    log_info "Static IP:   $static_ip"
    log_info "Access:      http://$static_ip"
    echo ""
    log_info "The VM is starting up (2-3 minutes for Docker + containers)."
    log_info "Monitor progress:"
    log_info "  gcloud compute ssh $vm_name --zone=$zone --command='tail -f /var/log/startup-script.log'"
    echo ""
    log_info "Lifecycle commands:"
    log_info "  Status:  ./status-vm.sh"
    log_info "  SSH:     ./ssh-vm.sh"
    log_info "  Logs:    ./logs-vm.sh"
    log_info "  Stop:    ./stop-vm.sh"
    log_info "  Start:   ./start-vm.sh"
    log_info "  Delete:  ./delete-vm.sh"
    echo ""
}

main "$@"
