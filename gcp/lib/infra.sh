#!/bin/bash
# Infrastructure library for v2.7 automated deployment
# Provides idempotent operations, logging, dry-run, cleanup, cost display, and VM creation
#
# Usage:
#   source gcp/lib/infra.sh
#   # Functions available: log_info, log_warn, log_error, execute,
#   # register_resource, cleanup_on_failure, create_static_ip,
#   # create_firewall_rules, create_persistent_disk, display_cost_estimate,
#   # create_vm

# Globals
DRY_RUN="${DRY_RUN:-false}"
CREATED_RESOURCES=()

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# ============================================================================
# Logging functions
# ============================================================================

log_info() {
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "${GREEN}[$timestamp] [INFO]${NC} $1"
}

log_warn() {
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "${YELLOW}[$timestamp] [WARN]${NC} $1"
}

log_error() {
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "${RED}[$timestamp] [ERROR]${NC} $1" >&2
}

# ============================================================================
# Dry-run wrapper
# ============================================================================

execute() {
    local description="$1"
    shift
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] $description"
        log_info "  Command: $*"
        return 0
    fi
    log_info "$description"
    "$@"
}

# ============================================================================
# Resource tracking and cleanup
# ============================================================================

register_resource() {
    # Format: "type:name:location"
    CREATED_RESOURCES+=("$1")
}

cleanup_on_failure() {
    local exit_code=$?
    if [[ $exit_code -ne 0 && ${#CREATED_RESOURCES[@]} -gt 0 ]]; then
        log_error "Deployment failed (exit code $exit_code). Cleaning up..."
        for resource in "${CREATED_RESOURCES[@]}"; do
            IFS=':' read -r type name location <<< "$resource"
            case "$type" in
                vm)
                    log_info "Deleting VM $name in $location..."
                    gcloud compute instances delete "$name" --zone="$location" --quiet 2>/dev/null || true
                    ;;
                # Never delete disks (data loss risk) or IPs (reusable)
            esac
        done
        log_info "Cleanup complete"
    fi
}

# ============================================================================
# Idempotent infrastructure functions
# ============================================================================

create_static_ip() {
    local ip_name="$1"
    local region="$2"

    # Check if IP exists (always run, not wrapped in execute)
    local existing_ip
    existing_ip=$(gcloud compute addresses describe "$ip_name" --region="$region" --format="value(address)" --quiet 2>/dev/null || true)

    if [[ -n "$existing_ip" ]]; then
        log_info "Static IP '$ip_name' already exists: $existing_ip"
        echo "$existing_ip"
        return 0
    fi

    # Create new IP
    execute "Creating static IP '$ip_name' in region $region" \
        gcloud compute addresses create "$ip_name" \
        --region="$region" \
        --quiet

    # Get the created IP (always run, not wrapped in execute)
    existing_ip=$(gcloud compute addresses describe "$ip_name" --region="$region" --format="value(address)" --quiet 2>/dev/null || true)
    log_info "Static IP created: $existing_ip"
    echo "$existing_ip"
}

create_firewall_rules() {
    local resource_prefix="$1"
    local vm_tag="${resource_prefix}-vm"

    # HTTP rule (port 80)
    local http_rule="${resource_prefix}-allow-http"
    if gcloud compute firewall-rules describe "$http_rule" --quiet &>/dev/null; then
        log_info "Firewall rule '$http_rule' already exists"
    else
        execute "Creating firewall rule '$http_rule' (TCP:80)" \
            gcloud compute firewall-rules create "$http_rule" \
            --direction=INGRESS \
            --priority=1000 \
            --network=default \
            --action=ALLOW \
            --rules=tcp:80 \
            --source-ranges=0.0.0.0/0 \
            --target-tags="$vm_tag" \
            --quiet
    fi

    # SSH rule (port 22)
    local ssh_rule="${resource_prefix}-allow-ssh"
    if gcloud compute firewall-rules describe "$ssh_rule" --quiet &>/dev/null; then
        log_info "Firewall rule '$ssh_rule' already exists"
    else
        execute "Creating firewall rule '$ssh_rule' (TCP:22)" \
            gcloud compute firewall-rules create "$ssh_rule" \
            --direction=INGRESS \
            --priority=1000 \
            --network=default \
            --action=ALLOW \
            --rules=tcp:22 \
            --source-ranges=0.0.0.0/0 \
            --target-tags="$vm_tag" \
            --quiet
    fi
}

create_persistent_disk() {
    local disk_name="$1"
    local zone="$2"
    local disk_size_gb="$3"

    # Check if disk exists (always run, not wrapped in execute)
    if gcloud compute disks describe "$disk_name" --zone="$zone" --quiet &>/dev/null; then
        local existing_size
        existing_size=$(gcloud compute disks describe "$disk_name" --zone="$zone" --format="value(sizeGb)" --quiet 2>/dev/null || true)
        log_info "Persistent disk '$disk_name' already exists (${existing_size}GB)"
        return 0
    fi

    # Create new disk
    execute "Creating persistent SSD disk '$disk_name' (${disk_size_gb}GB)" \
        gcloud compute disks create "$disk_name" \
        --zone="$zone" \
        --type=pd-ssd \
        --size="${disk_size_gb}GB" \
        --quiet
    log_info "Persistent disk '$disk_name' created"
}

# ============================================================================
# Cost estimation display
# ============================================================================

display_cost_estimate() {
    local machine_type="$1"
    local spot_price_hourly="$2"
    local disk_size_gb="$3"

    # Calculate monthly costs using bc
    local monthly_spot
    monthly_spot=$(echo "$spot_price_hourly * 730" | bc -l)

    local disk_cost
    disk_cost=$(echo "$disk_size_gb * 0.17" | bc -l)

    local total
    total=$(echo "$monthly_spot + $disk_cost" | bc -l)

    # Display itemized breakdown
    echo ""
    echo "Cost Estimate:"
    printf "  %-30s \$%.4f/hour  ~\$%.2f/month (730 hours)\n" "VM (Spot):" "$spot_price_hourly" "$monthly_spot"
    printf "  %-30s \$0.17/GB/month  ~\$%.2f/month (%d GB)\n" "Persistent disk (SSD):" "$disk_cost" "$disk_size_gb"
    printf "  %-30s %s\n" "Static IP (in use):" "Free"
    echo "  ──────────────────────────────────────────"
    printf "  %-30s ~\$%.2f/month\n" "Estimated total:" "$total"
    echo ""
    echo "Note: Spot pricing fluctuates. VM may be preempted with 30s notice."
    echo ""
}

# ============================================================================
# VM creation
# ============================================================================

create_vm() {
    local vm_name="$1"
    local zone="$2"
    local machine_type="$3"
    local static_ip="$4"
    local disk_name="$5"
    local startup_script_path="$6"
    local shutdown_script_path="$7"
    local resource_prefix="$8"

    # Check if VM exists (always run, not wrapped in execute)
    if gcloud compute instances describe "$vm_name" --zone="$zone" --quiet &>/dev/null; then
        log_error "VM '$vm_name' already exists in zone $zone"
        return 1
    fi

    # Create VM
    execute "Creating Spot VM '$vm_name' with machine type $machine_type" \
        gcloud compute instances create "$vm_name" \
        --zone="$zone" \
        --machine-type="$machine_type" \
        --provisioning-model=SPOT \
        --instance-termination-action=STOP \
        --tags="${resource_prefix}-vm" \
        --network-interface="address=$static_ip" \
        --boot-disk-size=20GB \
        --boot-disk-type=pd-balanced \
        --image-family=debian-12 \
        --image-project=debian-cloud \
        --disk="name=${disk_name},mode=rw,device-name=data-disk,boot=no,auto-delete=no" \
        --metadata-from-file=startup-script="${startup_script_path}",shutdown-script="${shutdown_script_path}" \
        --quiet

    # Register for cleanup on failure
    register_resource "vm:$vm_name:$zone"
    log_info "VM '$vm_name' created successfully"
}
