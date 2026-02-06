#!/bin/bash
# Pricing query library for v2.7 automated deployment
# Queries spot pricing and provides region selection functions
#
# Usage:
#   source gcp/lib/pricing.sh
#   region=$(get_cheapest_region 8 32)
#   zone=$(get_cheapest_zone 8 32)

get_cheapest_region() {
    local cpu_cores="${1:-8}"
    local ram_gb="${2:-32}"
    local script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

    local pricing_json
    pricing_json=$(python3 "$script_dir/query_pricing.py" \
        --cpu-cores "$cpu_cores" \
        --ram-gb "$ram_gb" 2>/dev/null) || {
        echo "ERROR: Pricing query failed" >&2
        return 1
    }

    # Return first (cheapest) region
    echo "$pricing_json" | jq -r '.[0].region'
}

get_cheapest_zone() {
    local cpu_cores="${1:-8}"
    local ram_gb="${2:-32}"
    local script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

    local pricing_json
    pricing_json=$(python3 "$script_dir/query_pricing.py" \
        --cpu-cores "$cpu_cores" \
        --ram-gb "$ram_gb" 2>/dev/null) || {
        echo "ERROR: Pricing query failed" >&2
        return 1
    }

    # Return first (cheapest) zone
    echo "$pricing_json" | jq -r '.[0].zone'
}

get_pricing_table() {
    local cpu_cores="${1:-8}"
    local ram_gb="${2:-32}"
    local script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

    # Return full JSON for display/processing
    python3 "$script_dir/query_pricing.py" \
        --cpu-cores "$cpu_cores" \
        --ram-gb "$ram_gb"
}
