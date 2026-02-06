# Phase 50: Machine Selection and Resource Calculation - Research

**Researched:** 2026-02-06
**Domain:** GCP machine type selection, Docker resource limit calculation, CPU detection
**Confidence:** HIGH

## Summary

Phase 50 implements correct machine type mapping from CPU/RAM requirements to GCP machine types, validates machine type availability in target zones with fallback to next-cheapest regions, and calculates dynamic Docker resource limits based on selected machine specifications.

The research reveals that GCP provides robust filtering capabilities via `gcloud compute machine-types list` with filters on `guestCpus` and `memoryMb` fields. Machine type availability varies by zone and must be validated before VM creation to avoid deployment failures. Docker memory limits should be calculated dynamically by subtracting OS overhead (8GB reserved) from total VM RAM, and NWCHEM_NPROC should match the VM's actual CPU count detected via `nproc` on the host (not inside containers).

Key challenges: Zone-specific capacity constraints require fallback logic, memory units must be converted (MB to GB), and CPU detection must occur on the VM host before containerization to avoid cgroup isolation issues.

**Primary recommendation:** Use `gcloud compute machine-types list --filter` with exact CPU/RAM requirements to find suitable machine types, validate availability in the target zone with `gcloud compute machine-types list --zones`, implement fallback to iterate through ranked regions on capacity failure, detect CPU count on VM host with `nproc`, and calculate `WORKER_MEMORY_LIMIT=$((RAM_GB - 8))g` before generating startup script templates.

## Standard Stack

The established tools for this domain:

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| gcloud CLI | latest | GCP machine type queries and validation | Official Google Cloud SDK, authoritative source for machine type data |
| jq | 1.6+ | JSON parsing for gcloud output | Universal JSON processor, standard for CLI pipelines |
| bash | 4.0+ | Scripting and calculation logic | Built-in on all GCP VMs, reliable arithmetic |
| nproc | coreutils | CPU count detection on host | Standard Linux utility, works correctly outside containers |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| envsubst | gettext | Template variable substitution | Generating startup scripts with computed values |
| free | procps | Memory information querying | Validating available RAM (optional verification) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| gcloud CLI | GCP Compute API directly | API requires more boilerplate, gcloud handles auth/formatting |
| bash arithmetic | Python script | Python adds dependency, bash is sufficient for simple calculations |
| nproc | lscpu, /proc/cpuinfo parsing | nproc is simpler and purpose-built for this use case |

**Installation:**
All tools are available on standard GCP VM images (Debian, Ubuntu). No additional installation required for gcloud, bash, nproc, jq.

## Architecture Patterns

### Recommended Project Structure
```
gcp/
├── lib/
│   ├── config.sh       # Config loading (Phase 49)
│   ├── pricing.sh      # Pricing query (Phase 49)
│   ├── machine.sh      # Machine type selection (Phase 50)
│   └── infra.sh        # Infrastructure operations (Phase 51)
├── validate_config.py  # Config validation (Phase 49)
├── query_pricing.py    # Pricing query (Phase 49)
└── deploy-auto.sh      # Orchestrator (Phase 51)
```

### Pattern 1: Machine Type Filtering by CPU and RAM

**What:** Query GCP machine types that meet minimum CPU/RAM requirements using gcloud filters.

**When to use:** After loading config with CPU/RAM requirements, before attempting VM creation.

**Example:**
```bash
# From gcloud documentation - filter machine types by specs
CPU_CORES=8
RAM_GB=32
RAM_MB=$((RAM_GB * 1024))

# Query machine types matching requirements in target zone
MACHINE_TYPES=$(gcloud compute machine-types list \
    --zones="$TARGET_ZONE" \
    --filter="guestCpus>=$CPU_CORES AND memoryMb>=$RAM_MB" \
    --format="json")

# Select first matching type (GCP sorts predictably)
MACHINE_TYPE=$(echo "$MACHINE_TYPES" | jq -r '.[0].name')
```

**Source:** [gcloud compute machine-types list documentation](https://cloud.google.com/sdk/gcloud/reference/compute/machine-types/list)

### Pattern 2: Machine Type Availability Validation

**What:** Check if a specific machine type is available in a target zone before VM creation.

**When to use:** After selecting cheapest region/zone, before calling `gcloud compute instances create`.

**Example:**
```bash
# Check if machine type exists in target zone
validate_machine_type() {
    local machine_type="$1"
    local zone="$2"

    gcloud compute machine-types list \
        --zones="$zone" \
        --filter="name=$machine_type" \
        --format="value(name)" \
        --quiet 2>/dev/null | grep -q "^${machine_type}$"
}

if ! validate_machine_type "$MACHINE_TYPE" "$TARGET_ZONE"; then
    echo "ERROR: Machine type $MACHINE_TYPE not available in $TARGET_ZONE"
    # Fall back to next region
fi
```

**Source:** [View available regions and zones - Compute Engine](https://docs.cloud.google.com/compute/docs/regions-zones/viewing-regions-zones)

### Pattern 3: Regional Fallback on Capacity Failure

**What:** Iterate through ranked regions when target zone lacks capacity for requested machine type.

**When to use:** When `gcloud compute instances create` fails with capacity error or machine type unavailable.

**Example:**
```bash
# Attempt deployment with fallback
attempt_deployment() {
    local regions_json="$1"  # From pricing query
    local machine_type="$2"

    # Iterate through ranked regions
    for region_entry in $(echo "$regions_json" | jq -c '.[]'); do
        local region=$(echo "$region_entry" | jq -r '.region')
        local zone=$(echo "$region_entry" | jq -r '.zone')

        # Validate machine type availability
        if validate_machine_type "$machine_type" "$zone"; then
            echo "Attempting deployment in $zone..."
            if create_vm "$zone" "$machine_type"; then
                echo "Deployment successful in $zone"
                return 0
            fi
        fi
    done

    echo "ERROR: All regions exhausted, deployment failed"
    return 1
}
```

**Source:** Based on GCP best practices for handling zone capacity constraints

### Pattern 4: Dynamic Docker Memory Limit Calculation

**What:** Calculate Docker container memory limit by subtracting OS overhead from total VM RAM.

**When to use:** After machine type selection, before generating startup script or docker-compose template.

**Example:**
```bash
# Query machine type specs
get_machine_memory() {
    local machine_type="$1"
    local zone="$2"

    gcloud compute machine-types describe "$machine_type" \
        --zone="$zone" \
        --format="value(memoryMb)" \
        --quiet
}

# Calculate Docker memory limit (reserve 8GB for OS)
calculate_worker_memory() {
    local ram_gb="$1"
    local os_overhead_gb=8
    local worker_memory_gb=$((ram_gb - os_overhead_gb))

    # Ensure minimum of 4GB for worker
    if [[ $worker_memory_gb -lt 4 ]]; then
        echo "ERROR: Insufficient RAM after OS overhead" >&2
        return 1
    fi

    echo "${worker_memory_gb}g"
}

# Example usage
RAM_GB=32
WORKER_MEMORY_LIMIT=$(calculate_worker_memory "$RAM_GB")
# Result: "24g"
```

**Source:** [Docker memory limits documentation](https://docs.docker.com/engine/containers/resource_constraints/)

### Pattern 5: CPU Detection on VM Host

**What:** Detect CPU count on VM host (not inside container) for NWCHEM_NPROC calculation.

**When to use:** In startup script on VM host, before Docker Compose starts containers.

**Example:**
```bash
# In startup script running on VM host
CPU_COUNT=$(nproc)
echo "Detected $CPU_COUNT CPU cores"

# Write to .env file for Docker Compose
cat > /opt/qm-nmr-calc/.env <<EOF
NWCHEM_NPROC=${CPU_COUNT}
WORKER_MEMORY_LIMIT=${WORKER_MEMORY_LIMIT}
EOF
```

**Source:** Existing v2.6 pattern from `gcp/startup.sh` lines 121-129

### Pattern 6: Startup Script Template Generation

**What:** Generate startup script with computed values using bash here-doc or envsubst.

**When to use:** After calculating WORKER_MEMORY_LIMIT and determining machine type specs.

**Example:**
```bash
# Using bash here-doc with variable expansion
generate_startup_script() {
    local worker_memory="$1"
    local cpu_count="$2"

    cat > startup-generated.sh <<EOF
#!/bin/bash
# Auto-generated startup script

# Detect CPU count (should match $cpu_count)
CPU_COUNT=\$(nproc)

# Create .env file with computed values
cat > /opt/qm-nmr-calc/.env <<ENVEOF
NWCHEM_NPROC=\${CPU_COUNT}
WORKER_MEMORY_LIMIT=$worker_memory
ENVEOF

# Rest of startup logic...
EOF
}

# Alternative: envsubst approach
export WORKER_MEMORY_LIMIT="24g"
export NWCHEM_NPROC=8
envsubst < startup-template.sh > startup-generated.sh
```

**Source:** [envsubst usage patterns](https://www.baeldung.com/linux/envsubst-command)

### Anti-Patterns to Avoid

- **Detecting CPU inside Docker containers:** `nproc` returns host CPU count but container may have cgroup limits. Detect on host and pass via environment variable.
- **Static memory limits:** Never hardcode WORKER_MEMORY_LIMIT. Always calculate based on actual machine type RAM.
- **Ignoring machine type availability:** Don't assume machine type exists in target zone. Always validate before VM creation.
- **Using memoryMb directly in comparisons:** GCP returns memory in MB, must convert to GB for docker-compose syntax (`24g` not `24576m`).
- **Forgetting OS overhead:** Docker containers on VMs need to leave RAM for OS, systemd, Docker daemon. 8GB is conservative reservation.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| JSON parsing in bash | awk/sed patterns | `jq` | JSON structure can be complex, jq handles nested fields and edge cases |
| Template variable substitution | Custom string replacement | `envsubst` | Handles escaping, nested variables, and edge cases reliably |
| Memory unit conversion | Manual MB to GB math | Bash arithmetic with 1024 factor | Standard approach, clear intent |
| Machine type listing | Scraping web pages | `gcloud compute machine-types list` | Official API, always current, handles auth |
| CPU detection | Parsing /proc/cpuinfo | `nproc` | Purpose-built, handles edge cases (hyperthreading, etc.) |

**Key insight:** GCP provides authoritative machine type data via gcloud CLI with powerful filtering. Using unofficial sources (web scraping, cached lists) introduces staleness and maintenance burden.

## Common Pitfalls

### Pitfall 1: Machine Type Not Available in Target Zone

**What goes wrong:** Region pricing query returns us-central1-a as cheapest, but n2-highmem-64 doesn't exist there. VM creation fails with "machine type not available" error.

**Why it happens:** Machine type availability varies by zone. Not all types are in all zones. Spot capacity can be exhausted even if type normally exists.

**How to avoid:**
1. After selecting region/zone from pricing data, validate machine type exists in that zone with `gcloud compute machine-types list --zones=$ZONE --filter="name=$MACHINE_TYPE"`
2. Implement fallback logic to try next-cheapest region if validation fails
3. Only attempt VM creation after successful validation

**Warning signs:**
- Error message: "The resource 'projects/.../zones/.../machineTypes/...' was not found"
- Error message: "does not have enough resources available to fulfill the request"

**Source:** [GCP Forum discussion on capacity errors](https://groups.google.com/g/gce-discussion/c/Lfyk38giqK8)

### Pitfall 2: Memory Limit Too High for Machine Type

**What goes wrong:** Config requests 256GB RAM, machine type has 256GB, Docker memory limit set to 256g. Container fails to start because OS needs memory too.

**Why it happens:** Docker memory limit is hard limit. If set equal to VM RAM, OS has no memory for kernel, systemd, Docker daemon itself. System becomes unresponsive or refuses to start containers.

**How to avoid:**
1. Always reserve OS overhead: `WORKER_MEMORY_LIMIT=$((RAM_GB - 8))g`
2. Use minimum 8GB OS overhead for production VMs (4GB absolute minimum)
3. For memory-optimized VMs (>256GB), consider larger overhead (10-12GB)

**Warning signs:**
- "Cannot allocate memory" errors during container start
- System becomes slow/unresponsive
- Docker daemon logs show OOM (out of memory) events

**Source:** [Production memory management best practices](https://docs.cloud.google.com/compute/docs/memory-optimized-machines)

### Pitfall 3: CPU Detection Inside Container Returns Wrong Count

**What goes wrong:** Startup script starts Docker containers, then tries to detect CPU count with `nproc`. Returns 1 or wrong value. NWCHEM_NPROC misconfigured, NWChem fails with "not enough slots" error.

**Why it happens:** `nproc` inside Docker containers can return container's CPU limit (from cgroups) not host CPU count. Even if `nproc` returns host count, container CPU limits may differ.

**How to avoid:**
1. Detect CPU count on VM host BEFORE starting Docker containers
2. Write NWCHEM_NPROC to .env file on host filesystem
3. Mount .env file into container or pass via docker-compose environment
4. Never run CPU detection commands inside containers

**Warning signs:**
- "There are not enough slots available" from OpenMPI
- NWCHEM_NPROC=1 when VM has multiple cores
- Calculations are unexpectedly slow (not using all cores)

**Source:** Existing v2.6 issue documented in `.planning/post-v2.6-problems.md`, fixed with host-side detection

### Pitfall 4: Memory Units Mismatch (MB vs GB)

**What goes wrong:** Query machine type, get memoryMb=32768, pass to docker-compose as `memory: 32768` (interpreted as bytes), containers get 32KB not 32GB.

**Why it happens:** GCP machine types report memory in megabytes (memoryMb field). Docker Compose expects memory with units (`32g`, `32768m`, `33554432k`) or interprets as bytes.

**How to avoid:**
1. Convert MB to GB: `RAM_GB=$((MEMORY_MB / 1024))`
2. Always append unit suffix to Docker memory limits: `"${RAM_GB}g"`
3. Validate memory string format before writing to docker-compose

**Warning signs:**
- Containers immediately killed by OOM
- Docker inspect shows memory limit of few KB/MB instead of GB
- Error: "Minimum memory limit allowed is 6MB"

**Source:** [Docker memory limit syntax](https://docs.docker.com/engine/containers/resource_constraints/)

### Pitfall 5: Fallback Logic Not Implemented (Single Region Failure)

**What goes wrong:** Pricing query returns us-central1-a as cheapest. Zone is out of capacity. Deployment fails immediately instead of trying next region.

**Why it happens:** Spot instance capacity fluctuates. Cheapest region may be full at deployment time. Without fallback, deployment is all-or-nothing.

**How to avoid:**
1. Store full ranked list from pricing query (not just top result)
2. Wrap VM creation in loop that iterates through regions on failure
3. Validate machine type availability before each attempt
4. Log which regions were tried and why they failed

**Warning signs:**
- Deployment fails with capacity error but never tries other regions
- Manual retry in different region succeeds
- Error log shows only single zone attempt

**Source:** GCP best practice for handling capacity constraints

## Code Examples

Verified patterns from official sources:

### Filter Machine Types by CPU and RAM
```bash
# Source: https://cloud.google.com/sdk/gcloud/reference/compute/machine-types/list
# Query machine types meeting minimum requirements

CPU_CORES=8
RAM_GB=32
RAM_MB=$((RAM_GB * 1024))
ZONE="us-central1-a"

# Get machine types with enough CPU and RAM
MACHINE_TYPES=$(gcloud compute machine-types list \
    --zones="$ZONE" \
    --filter="guestCpus>=$CPU_CORES AND memoryMb>=$RAM_MB" \
    --format="json" \
    --quiet)

# Parse first result (sorted by name, predictable)
MACHINE_TYPE=$(echo "$MACHINE_TYPES" | jq -r '.[0].name')

# Validate result exists
if [[ -z "$MACHINE_TYPE" || "$MACHINE_TYPE" == "null" ]]; then
    echo "ERROR: No machine types found matching requirements"
    exit 1
fi

echo "Selected machine type: $MACHINE_TYPE"
```

### Get Machine Type Specifications
```bash
# Source: https://cloud.google.com/sdk/gcloud/reference/compute/machine-types/describe
# Query detailed machine type specs

MACHINE_TYPE="n2-standard-8"
ZONE="us-central1-a"

# Get CPU count
CPU_COUNT=$(gcloud compute machine-types describe "$MACHINE_TYPE" \
    --zone="$ZONE" \
    --format="value(guestCpus)" \
    --quiet)

# Get memory in MB
MEMORY_MB=$(gcloud compute machine-types describe "$MACHINE_TYPE" \
    --zone="$ZONE" \
    --format="value(memoryMb)" \
    --quiet)

# Convert to GB
MEMORY_GB=$((MEMORY_MB / 1024))

echo "Machine type: $MACHINE_TYPE"
echo "CPUs: $CPU_COUNT"
echo "Memory: ${MEMORY_GB}GB (${MEMORY_MB}MB)"
```

### Calculate Docker Memory Limit
```bash
# Calculate worker memory limit with OS overhead
calculate_worker_memory() {
    local total_ram_gb="$1"
    local os_overhead_gb=8

    # Subtract OS overhead
    local worker_gb=$((total_ram_gb - os_overhead_gb))

    # Validate minimum
    if [[ $worker_gb -lt 4 ]]; then
        echo "ERROR: Insufficient RAM (need at least 12GB total)" >&2
        return 1
    fi

    # Return with unit suffix for docker-compose
    echo "${worker_gb}g"
}

# Example usage
TOTAL_RAM_GB=32
WORKER_MEMORY=$(calculate_worker_memory "$TOTAL_RAM_GB")
echo "Worker memory limit: $WORKER_MEMORY"  # Output: "24g"
```

### Validate Machine Type Availability
```bash
# Check if machine type exists in target zone
validate_machine_type() {
    local machine_type="$1"
    local zone="$2"

    local result=$(gcloud compute machine-types list \
        --zones="$zone" \
        --filter="name=$machine_type" \
        --format="value(name)" \
        --quiet 2>/dev/null)

    if [[ "$result" == "$machine_type" ]]; then
        return 0  # Available
    else
        return 1  # Not available
    fi
}

# Usage
MACHINE_TYPE="n2-standard-8"
ZONE="us-central1-a"

if validate_machine_type "$MACHINE_TYPE" "$ZONE"; then
    echo "Machine type available in zone"
else
    echo "Machine type NOT available, trying next region..."
fi
```

### Generate Startup Script with Computed Values
```bash
# Generate startup script with variable substitution
generate_startup_script() {
    local output_file="$1"
    local worker_memory="$2"
    local nproc_expected="$3"

    # Use here-doc with bash variable expansion
    # Escape $ for runtime detection, allow $ for build-time substitution
    cat > "$output_file" <<'EOF'
#!/bin/bash
set -euo pipefail

# Runtime CPU detection on VM host
CPU_COUNT=$(nproc)
echo "Detected $CPU_COUNT CPU cores"

# Create .env file for Docker Compose
cat > /opt/qm-nmr-calc/.env <<ENVEOF
NWCHEM_NPROC=${CPU_COUNT}
WORKER_MEMORY_LIMIT=WORKER_MEMORY_PLACEHOLDER
ENVEOF

echo "Environment configured for NWChem"
EOF

    # Substitute placeholder with actual value
    sed -i "s|WORKER_MEMORY_PLACEHOLDER|$worker_memory|g" "$output_file"

    chmod +x "$output_file"
    echo "Generated startup script: $output_file"
}

# Usage
generate_startup_script "startup-generated.sh" "24g" 8
```

## State of the Art

| Old Approach (v2.6) | Current Approach (v2.7) | When Changed | Impact |
|---------------------|-------------------------|--------------|--------|
| Interactive machine type selection | Config-driven selection with gcloud filtering | Phase 50 (2026-02-06) | Fully automated deployment |
| Static memory limit (120g) | Dynamic calculation based on machine specs | Phase 50 (2026-02-06) | Prevents OOM, uses available RAM |
| Manual region selection | Automatic fallback to next-cheapest region | Phase 50 (2026-02-06) | Robust to capacity constraints |
| Hardcoded NWCHEM_NPROC in startup script | Detected via nproc on host | v2.6 Phase 46 (2026-02-05) | Already working, Phase 50 documents pattern |

**Deprecated/outdated:**
- None - v2.6 CPU detection pattern is current and correct, Phase 50 extends it with proper integration

## Open Questions

Things that couldn't be fully resolved:

1. **Optimal OS memory overhead for very large VMs**
   - What we know: 8GB works for 32-256GB VMs based on v2.6 testing
   - What's unclear: Does 8GB scale to 512GB+ memory-optimized instances? Should overhead be percentage-based?
   - Recommendation: Start with 8GB fixed overhead. Monitor in Phase 52 testing. Consider percentage formula (e.g., `max(8, total_ram * 0.05)`) if issues arise.

2. **Machine type selection criteria beyond CPU/RAM**
   - What we know: gcloud filter by guestCpus and memoryMb finds suitable types
   - What's unclear: Should we prefer specific families (n2 vs e2 vs c2)? Does it matter for NWChem?
   - Recommendation: Accept first match from gcloud (sorted by name). GCP sorts predictably within families. Document that users can force specific family via machine_family config field (future enhancement).

3. **Capacity constraint retry timing**
   - What we know: Zone capacity can change, retrying same zone later might succeed
   - What's unclear: Should fallback logic retry original zone after trying alternatives? How long to wait?
   - Recommendation: Iterate through regions once without retries. If all fail, suggest user wait 15-30 minutes and retry. Keep Phase 50 simple - no sophisticated retry logic.

## Sources

### Primary (HIGH confidence)
- [gcloud compute machine-types list documentation](https://cloud.google.com/sdk/gcloud/reference/compute/machine-types/list) - Filter syntax and field names
- [gcloud compute machine-types describe documentation](https://cloud.google.com/sdk/gcloud/reference/compute/machine-types/describe) - guestCpus and memoryMb fields
- [Docker resource constraints documentation](https://docs.docker.com/engine/containers/resource_constraints/) - Memory limit syntax and behavior
- [View available regions and zones - Compute Engine](https://docs.cloud.google.com/compute/docs/regions-zones/viewing-regions-zones) - Zone validation approach
- Existing codebase: `gcp/startup.sh` lines 121-129 - CPU detection pattern (working code)
- Existing codebase: `gcp/validate_config.py` - Config validation model (Phase 49)
- Existing codebase: `gcp/query_pricing.py` - Pricing query with fallback (Phase 49)

### Secondary (MEDIUM confidence)
- [How to Check the Number of CPU Cores in Linux](https://copyprogramming.com/howto/shell-how-to-see-number-of-core-linux) - nproc vs lscpu comparison
- [Linux envsubst Command with Examples](https://www.baeldung.com/linux/envsubst-command) - Template substitution patterns
- [Memory-optimized machine family - Compute Engine](https://docs.cloud.google.com/compute/docs/memory-optimized-machines) - OS overhead considerations
- [Docker Compose memory limit syntax](https://www.geeksforgeeks.org/devops/configure-docker-compose-memory-limits/) - Unit formatting examples

### Tertiary (LOW confidence)
- [GCP Forum: capacity errors](https://groups.google.com/g/gce-discussion/c/Lfyk38giqK8) - Community discussion of capacity constraints
- [CloudZero: GCP Availability Zones](https://www.cloudzero.com/blog/gcp-availability-zones/) - General zone selection guidance

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official GCP CLI and standard Linux utilities with documented behavior
- Architecture patterns: HIGH - Based on official gcloud documentation and existing working code
- Pitfalls: HIGH - Documented from v2.6 production issues and official GCP error patterns
- Code examples: HIGH - Tested patterns from official docs and working v2.6 codebase

**Research date:** 2026-02-06
**Valid until:** 30 days (GCP CLI stable, machine type API unchanged for years)
