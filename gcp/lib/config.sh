#!/bin/bash
# Config loading library for v2.7 automated deployment
# Sources validated TOML config as bash environment variables
#
# Usage:
#   source gcp/lib/config.sh
#   load_config [path/to/config.toml]

load_config() {
    local config_path="${1:-gcp/config.toml}"
    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

    if [[ ! -f "$config_path" ]]; then
        echo "ERROR: Config file not found: $config_path" >&2
        echo "Create from template: cp gcp/config.toml.example gcp/config.toml" >&2
        return 1
    fi

    local exports
    exports=$(python3 "$script_dir/validate_config.py" --config "$config_path" 2>&1) || {
        echo "$exports" >&2
        return 1
    }

    eval "$exports"
}
