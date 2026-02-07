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

    # ── Step 3: Get ranked zones for machine selection ──
    log_info "Step 3/6: Finding machine type and ranked zones..."
    local machine_json
    machine_json=$(select_machine "$CPU_CORES" "$RAM_GB") || { log_error "Machine selection failed"; exit 1; }

    local machine_type first_zone first_region
    machine_type=$(echo "$machine_json" | jq -r '.machine_type')
    first_zone=$(echo "$machine_json" | jq -r '.zone')
    first_region=$(echo "$machine_json" | jq -r '.region')

    log_info "  Machine type: $machine_type"
    log_info "  Primary zone: $first_zone"
    echo ""

    # Get spot pricing for cost estimate
    local pricing_json spot_price_hourly
    pricing_json=$(get_pricing_table "$CPU_CORES" "$RAM_GB") || true
    spot_price_hourly=$(echo "$pricing_json" | jq -r '.[0].spot_price_hourly // "0.05"')

    # ── Step 4: Display cost estimate ──
    log_info "Step 4/6: Cost estimate..."
    display_cost_estimate "$machine_type" "$spot_price_hourly" "$DISK_SIZE_GB"
    echo ""

    # ── Step 5: Create shared infrastructure ──
    log_info "Step 5/6: Creating infrastructure..."
    local ip_name="${RESOURCE_PREFIX}-ip"
    local disk_name="${RESOURCE_PREFIX}-data"
    local vm_name="${RESOURCE_PREFIX}-vm"

    create_firewall_rules "$RESOURCE_PREFIX"
    echo ""

    # ── Step 6: Try zones until Spot VM succeeds ──
    log_info "Step 6/6: Creating Spot VM..."

    # Build zone list from ranked regions (primary first, then alternatives)
    local zones=()
    local zone_entry
    while IFS= read -r zone_entry; do
        [[ -n "$zone_entry" ]] && zones+=("$zone_entry")
    done < <(echo "$pricing_json" | jq -r '.[].zone' 2>/dev/null)

    # Ensure primary zone is first
    if [[ ${#zones[@]} -eq 0 ]]; then
        zones=("$first_zone")
    fi

    # Use existing shutdown.sh from gcp/
    local shutdown_script="$SCRIPT_DIR/shutdown.sh"
    if [[ ! -f "$shutdown_script" ]]; then
        log_error "Shutdown script not found: $shutdown_script"
        exit 1
    fi

    local zone region static_ip vm_created=false
    local prev_region=""
    for zone in "${zones[@]}"; do
        region="${zone%-*}"

        # Check machine type exists in this zone
        if [[ "$DRY_RUN" != "true" ]]; then
            if ! gcloud compute machine-types describe "$machine_type" --zone="$zone" --quiet &>/dev/null; then
                log_warn "Machine type $machine_type not available in $zone, skipping..."
                continue
            fi
        fi

        log_info "Trying zone: $zone (region: $region)"

        # Handle region change: static IPs are regional
        if [[ -n "$prev_region" && "$region" != "$prev_region" ]]; then
            log_info "Region changed ($prev_region → $region), migrating static IP..."
            delete_static_ip "$ip_name" "$prev_region"
        fi

        # Create regional IP and zonal disk
        static_ip=$(create_static_ip "$ip_name" "$region")
        create_persistent_disk "$disk_name" "$zone" "$DISK_SIZE_GB"

        # Generate dynamic startup script
        local startup_tmp
        startup_tmp=$(mktemp /tmp/qm-nmr-startup.XXXXXX)
        generate_startup "$CPU_CORES" "$RAM_GB" "$RESOURCE_PREFIX" "$DISK_SIZE_GB" > "$startup_tmp"

        # Try VM creation with error classification
        local vm_rc=0
        try_create_vm "$vm_name" "$zone" "$machine_type" "$static_ip" "$disk_name" \
            "$startup_tmp" "$shutdown_script" "$RESOURCE_PREFIX" || vm_rc=$?
        rm -f "$startup_tmp"

        case $vm_rc in
            0)
                vm_created=true
                break
                ;;
            1)
                # Retryable — clean up zonal disk and continue
                log_warn "Spot capacity exhausted in $zone, trying next zone..."
                if gcloud compute disks describe "$disk_name" --zone="$zone" --quiet &>/dev/null; then
                    log_info "Cleaning up disk in $zone for zone retry..."
                    gcloud compute disks delete "$disk_name" --zone="$zone" --quiet 2>/dev/null || true
                fi
                prev_region="$region"
                continue
                ;;
            2)
                # Fatal — stop retrying
                log_error "Fatal error creating VM — aborting zone retry"
                exit 1
                ;;
        esac
    done

    if [[ "$vm_created" != "true" ]]; then
        log_error "All zones exhausted — no Spot capacity for $machine_type"
        log_error "Try again later or use a different machine spec"
        exit 1
    fi

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
