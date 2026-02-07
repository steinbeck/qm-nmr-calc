#!/bin/bash
# Machine selection library for v2.7 automated deployment
# Wraps gcp/select_machine.py for bash consumption
#
# Usage:
#   source gcp/lib/machine.sh
#   result_json=$(select_machine 8 32)
#   eval $(get_docker_resources 8 32)
#   generate_startup 8 32 qm-nmr-calc 100 > startup.sh
#   get_machine_info 8 32

# Module-level path setup (computed once when sourced)
_MACHINE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
_MACHINE_PROJECT_ROOT="$(cd "$_MACHINE_SCRIPT_DIR/.." && pwd)"
export PYTHONPATH="$_MACHINE_PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

select_machine() {
    local cpu_cores="${1:-8}"
    local ram_gb="${2:-32}"

    local machine_json stderr_tmp
    stderr_tmp=$(mktemp)
    machine_json=$(python3 "$_MACHINE_SCRIPT_DIR/select_machine.py" \
        --cpu-cores "$cpu_cores" \
        --ram-gb "$ram_gb" 2>"$stderr_tmp") || {
        echo "ERROR: Machine selection failed: $(cat "$stderr_tmp")" >&2
        rm -f "$stderr_tmp"
        return 1
    }
    # Show warnings to user but keep them out of JSON
    [[ -s "$stderr_tmp" ]] && cat "$stderr_tmp" >&2
    rm -f "$stderr_tmp"

    echo "$machine_json"
}

get_docker_resources() {
    local cpu_cores="${1:-8}"
    local ram_gb="${2:-32}"

    local machine_json
    machine_json=$(python3 "$_MACHINE_SCRIPT_DIR/select_machine.py" \
        --cpu-cores "$cpu_cores" \
        --ram-gb "$ram_gb" 2>/dev/null) || {
        echo "ERROR: Machine selection failed" >&2
        return 1
    }

    # Extract worker_memory_limit and nwchem_nproc for eval
    local worker_memory_limit
    local nwchem_nproc

    worker_memory_limit=$(echo "$machine_json" | jq -r '.resources.worker_memory_limit')
    nwchem_nproc=$(echo "$machine_json" | jq -r '.resources.nwchem_nproc')

    # Output for eval pattern
    echo "WORKER_MEMORY_LIMIT=$worker_memory_limit"
    echo "NWCHEM_NPROC=$nwchem_nproc"
}

generate_startup() {
    local cpu_cores="${1:-8}"
    local ram_gb="${2:-32}"
    local resource_prefix="${3:-qm-nmr-calc}"
    local disk_size_gb="${4:-100}"

    # Call Python CLI with --generate-startup-script flag
    python3 "$_MACHINE_SCRIPT_DIR/select_machine.py" \
        --cpu-cores "$cpu_cores" \
        --ram-gb "$ram_gb" \
        --resource-prefix "$resource_prefix" \
        --disk-size-gb "$disk_size_gb" \
        --generate-startup-script
}

get_machine_info() {
    local cpu_cores="${1:-8}"
    local ram_gb="${2:-32}"

    local machine_json
    machine_json=$(python3 "$_MACHINE_SCRIPT_DIR/select_machine.py" \
        --cpu-cores "$cpu_cores" \
        --ram-gb "$ram_gb" 2>/dev/null) || {
        echo "ERROR: Machine selection failed" >&2
        return 1
    }

    # Extract fields for human-readable output
    local machine_type
    local zone
    local total_ram_gb
    local total_cpus

    machine_type=$(echo "$machine_json" | jq -r '.machine_type')
    zone=$(echo "$machine_json" | jq -r '.zone')
    total_ram_gb=$(echo "$machine_json" | jq -r '.resources.total_ram_gb')
    total_cpus=$(echo "$machine_json" | jq -r '.resources.total_cpus')

    # Human-readable format
    echo "Machine type: $machine_type in $zone (${total_ram_gb}GB RAM, ${total_cpus} CPUs)"
}
