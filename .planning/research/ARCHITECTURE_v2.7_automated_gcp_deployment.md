# Architecture Patterns: Config-Driven GCP Spot Deployment

**Domain:** Non-interactive cloud deployment automation
**Researched:** 2026-02-06
**Confidence:** MEDIUM (GCP pricing API limitations, third-party tooling required)

## Executive Summary

This architecture research addresses how to replace interactive GCP deployment scripts with a fully automated, config-driven system that finds the cheapest spot instance across all regions and deploys without prompts.

**Key findings:**
1. **GCP Billing API is not practical for real-time pricing** - No direct SKU lookup by machine type, requires filtering hundreds of SKUs by name matching
2. **Third-party pricing data is necessary** - CloudPrice.net and similar services aggregate spot pricing, but no official API
3. **Config-driven deployment requires early validation** - Parse config → validate specs → find cheapest region → deploy in single pass
4. **Memory/CPU automation is straightforward** - `nproc` works on host, write to `.env` file that Docker Compose reads
5. **Script organization: Modular functions, single entry point** - Break into sourced libraries for testability, single `deploy-auto.sh` orchestrates

## Recommended Architecture

### High-Level Flow

```
┌─────────────────┐
│ deploy-auto.sh  │  Single entry point, non-interactive
└────────┬────────┘
         │
         ├─> Load config.yaml (CPU cores, RAM, GCP project)
         │
         ├─> Validate config (CPU/RAM feasible? Project exists?)
         │
         ├─> Query pricing data (find cheapest region for specs)
         │
         ├─> Calculate Docker memory limits and NWCHEM_NPROC
         │
         ├─> Generate deployment parameters
         │
         ├─> Create infrastructure (static IP, firewall, disk) if missing
         │
         ├─> Deploy VM with startup script
         │
         └─> Output deployment summary (IP, cost estimate, SSH command)
```

### Component Boundaries

| Component | Responsibility | Dependencies |
|-----------|---------------|--------------|
| `deploy-auto.sh` | Orchestration, error handling, output | All libraries |
| `lib/config.sh` | Parse YAML config, validate specs | `yq` or built-in parser |
| `lib/pricing.sh` | Query pricing data, find cheapest region | Third-party pricing API or scraping |
| `lib/machine.sh` | Map specs to machine types, calculate limits | `gcloud compute machine-types list` |
| `lib/infra.sh` | Create/verify infrastructure resources | Existing `setup-infrastructure.sh` logic |
| `lib/deploy.sh` | VM creation with computed parameters | Existing `deploy-vm.sh` logic |
| `startup-template.sh` | Template for VM startup script | `.env` generation logic |

### Data Flow

1. **Config → Specs**: User provides `cpu_cores: 64` and `ram_gb: 256` in `config.yaml`
2. **Specs → Machine Types**: Query `gcloud compute machine-types list --filter="guestCpus=64 AND memoryMb>=262144"` across all regions
3. **Machine Types → Pricing**: Match machine types to pricing data (third-party API or cached data)
4. **Pricing → Region Selection**: Sort by price, select cheapest with availability
5. **Specs → Docker Limits**: Calculate `WORKER_MEMORY_LIMIT=$((ram_gb - 8))g` (leave 8GB for OS), `NWCHEM_NPROC=$cpu_cores`
6. **All Parameters → Deployment**: Pass to infrastructure creation and VM deployment

## Integration with Existing Scripts

### Refactor vs. New Script

**Recommendation: Extract libraries from existing scripts, create new orchestrator**

Rationale:
- Existing scripts have well-tested logic for infrastructure creation (`setup-infrastructure.sh`) and VM deployment (`deploy-vm.sh`)
- Interactive prompts are scattered throughout, hard to cleanly remove
- Modular approach allows reuse: `deploy-auto.sh` sources the same libraries as interactive scripts

### Reusable Components from Existing Scripts

| Existing Script | Extract to Library | Reuse in |
|-----------------|-------------------|----------|
| `setup-infrastructure.sh` | `lib/infra.sh` | Both interactive and automated |
| `deploy-vm.sh` | `lib/deploy.sh` | Both interactive and automated |
| `config.sh` | Keep as-is | All scripts |
| `startup.sh` | Template with variables | Generate dynamically |
| `docker-compose.gcp.yml` | Template with variables | Generate dynamically |

### Migration Path

1. **Phase 1: Extract functions** - Move infrastructure and deployment logic to `lib/` directory, keep existing scripts working
2. **Phase 2: Create orchestrator** - Write `deploy-auto.sh` that sources libraries and calls functions non-interactively
3. **Phase 3: Add config parsing** - Implement YAML config loader in `lib/config.sh`
4. **Phase 4: Integrate pricing** - Add `lib/pricing.sh` with third-party pricing data
5. **Phase 5: Deprecate interactive** - Document that `deploy-vm.sh` is now legacy, recommend `deploy-auto.sh`

## Config File Design

### Format: YAML (Not TOML or env)

**Rationale:**
- YAML is readable for non-developers (scientists deploying qm-nmr-calc)
- `yq` is a lightweight, single-binary parser available on GitHub releases
- `.env` files are not suitable (no nested structure, used for runtime config)
- TOML is overkill (no nested structure needed)
- JSON is not human-friendly for editing

### Schema

```yaml
# config.yaml - GCP automated deployment configuration
gcp:
  project_id: "my-gcp-project"
  resource_prefix: "qm-nmr-calc"

compute:
  # Desired CPU cores and RAM - script finds cheapest machine type
  cpu_cores: 64
  ram_gb: 256

  # Persistent disk for job data
  disk_size_gb: 100

  # Optional: restrict to specific regions (default: all regions)
  # regions: ["us-central1", "europe-west1"]

  # Optional: exclude regions (e.g., regulatory compliance)
  # exclude_regions: ["asia-*"]

deployment:
  # HTTP only (no domain, no HTTPS)
  http_only: true

  # Optional: custom image tags (default: latest)
  # image_tag: "v2.7.0"
```

### Validation Rules

```bash
# Required fields
- gcp.project_id (non-empty, no spaces)
- compute.cpu_cores (integer, 2-224)
- compute.ram_gb (integer, 4-3844)

# Constraints
- RAM must be >= 2GB per core (NWChem requirement)
- CPU cores must be available in at least one GCP machine type
- Disk size >= 10GB

# Warnings (non-fatal)
- CPU > 64: "High core counts may be expensive"
- RAM > 256GB: "High memory VMs have limited regional availability"
```

## GCP Pricing Query Architecture

### Problem: Official API is Impractical

**GCP Cloud Billing API limitations:**
1. No direct "get price for machine type X in region Y" endpoint
2. Requires filtering 500+ SKUs by display name string matching (e.g., "Spot Preemptible n2-standard-64 Instance Core running in Iowa")
3. Pricing is per-resource (CPU hours, memory hours) not per-instance
4. Requires computing total cost from multiple SKUs (CPU + RAM + disk + network)
5. Spot pricing is "dynamic" - API returns max price, not current market price

**Source:** [Get Google Cloud pricing information API](https://docs.cloud.google.com/billing/docs/how-to/get-pricing-information-api)

### Solution: Third-Party Pricing Aggregation

**Recommended approach: Use third-party pricing data**

Options ranked by reliability:

#### Option 1: CloudPrice.net JSON scraping (MEDIUM confidence)
- **URL:** `https://cloudprice.net/gcp/compute` provides real-time spot pricing
- **Pros:** Updated daily, covers all regions, includes spot pricing
- **Cons:** No official API, requires HTML/JSON parsing, may break on UI changes
- **Implementation:**
  ```bash
  # Fetch pricing data, cache for 24 hours
  curl -s 'https://cloudprice.net/gcp/compute' | \
    jq '.[] | select(.type == "spot") | {machine: .name, region: .region, price: .price_hour}'
  ```

**Sources:**
- [GCP Compute Engine Pricing - CloudPrice](https://cloudprice.net/gcp/compute)
- [GCP Region Pricing Comparison](https://cloudprice.net/gcp/regions)

#### Option 2: gcloud-pricing-calculator (LOW confidence)
- **GitHub:** [sagemathinc/gcloud-pricing-calculator](https://github.com/sagemathinc/gcloud-pricing-calculator)
- **Pros:** Open source, JavaScript library
- **Cons:** Not actively maintained, may have stale data, requires Node.js
- **Not recommended** due to maintenance uncertainty

#### Option 3: Manual price table (LOW confidence)
- **Approach:** Hard-code prices for common machine types in script
- **Pros:** No external dependencies
- **Cons:** Requires manual updates, quickly becomes stale
- **Only use as fallback** if pricing query fails

### Pricing Query Implementation Pattern

```bash
# lib/pricing.sh

get_cheapest_region() {
    local cpu_cores=$1
    local ram_gb=$2

    # Find machine types matching specs across all regions
    local machine_types=$(gcloud compute machine-types list \
        --filter="guestCpus=${cpu_cores} AND memoryMb>=$((ram_gb * 1024))" \
        --format="csv[no-heading](name,zone,guestCpus,memoryMb)")

    # Fetch pricing data (cached for 24 hours)
    local pricing_data=$(fetch_pricing_data)

    # Match machine types to pricing, find cheapest
    local cheapest_region=""
    local cheapest_price="999999"

    while IFS=',' read -r machine zone cpus mem; do
        local region=${zone%-*}
        local price=$(echo "$pricing_data" | jq -r \
            ".[] | select(.machine == \"$machine\" and .region == \"$region\") | .price")

        if (( $(echo "$price < $cheapest_price" | bc -l) )); then
            cheapest_price=$price
            cheapest_region=$region
            cheapest_zone=$zone
            cheapest_machine=$machine
        fi
    done <<< "$machine_types"

    echo "$cheapest_region $cheapest_zone $cheapest_machine $cheapest_price"
}
```

**Confidence level:** MEDIUM - Third-party pricing data is reliable but not guaranteed to be real-time or accurate.

## Docker Memory and CPU Configuration

### Automatic Limit Calculation

**Problem from v2.6:** Docker `nproc` returns 1 inside containers, memory limits were hardcoded to 8GB (insufficient for 64-core VM with 256GB RAM).

**Solution:** Calculate limits on VM host, write to `.env` file before starting containers.

### Implementation in Startup Script

```bash
#!/bin/bash
# startup-auto.sh (generated from template)

# Detect VM specifications
CPU_COUNT=$(nproc)  # Works correctly on host (not in container)
TOTAL_RAM_GB=$(free -g | awk '/^Mem:/ {print $2}')

# Calculate Docker limits
# Rule: Leave 8GB for OS, rest for worker container
WORKER_MEMORY_GB=$((TOTAL_RAM_GB - 8))
NWCHEM_NPROC=$CPU_COUNT

# Create .env file with computed values
cat > /opt/qm-nmr-calc/.env <<EOF
DOMAIN=
NWCHEM_NPROC=${NWCHEM_NPROC}
WORKER_MEMORY_LIMIT=${WORKER_MEMORY_GB}g
EOF

# docker-compose.gcp.yml template with variable substitution
cat > /opt/qm-nmr-calc/docker-compose.gcp.yml <<EOF
services:
  worker:
    volumes:
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info
      - NWCHEM_NPROC=${NWCHEM_NPROC}
    deploy:
      resources:
        limits:
          memory: ${WORKER_MEMORY_GB}g
EOF

# Start services
docker compose -f docker-compose.yml -f docker-compose.gcp.yml up -d
```

### Why This Works

1. **Host vs. Container CPU detection:** `nproc` on the VM host correctly returns 64 cores, but inside a container returns 1 (cgroups isolation issue)
2. **Memory calculation:** VM has 256GB total, leave 8GB for OS/overhead, allocate 248GB to worker container
3. **MPI oversubscribe:** Already fixed in v2.6 with `--oversubscribe` flag, still works with correct NWCHEM_NPROC
4. **No hardcoded limits:** All limits computed dynamically based on actual VM specs

**Source:** [Docker memory limits](https://docs.docker.com/engine/containers/resource_constraints/)

## Patterns to Follow

### Pattern 1: Config-First Validation

**What:** Validate configuration before making any GCP API calls

**When:** User runs `deploy-auto.sh`, before infrastructure creation

**Example:**
```bash
# lib/config.sh

validate_config() {
    local config_file=$1

    # Parse YAML
    local cpu_cores=$(yq e '.compute.cpu_cores' "$config_file")
    local ram_gb=$(yq e '.compute.ram_gb' "$config_file")
    local project_id=$(yq e '.gcp.project_id' "$config_file")

    # Validate CPU cores
    if [[ $cpu_cores -lt 2 || $cpu_cores -gt 224 ]]; then
        echo_error "CPU cores must be between 2 and 224 (got: $cpu_cores)"
        return 1
    fi

    # Validate RAM (must be >= 2GB per core for NWChem)
    local min_ram=$((cpu_cores * 2))
    if [[ $ram_gb -lt $min_ram ]]; then
        echo_error "RAM must be at least ${min_ram}GB for ${cpu_cores} cores (got: ${ram_gb}GB)"
        return 1
    fi

    # Validate GCP project exists
    if ! gcloud projects describe "$project_id" &>/dev/null; then
        echo_error "GCP project '$project_id' not found or no access"
        return 1
    fi

    echo_info "Configuration valid"
    return 0
}
```

### Pattern 2: Idempotent Infrastructure

**What:** Infrastructure creation should be safe to run multiple times (create if missing, skip if exists)

**When:** Every deployment, before VM creation

**Example:**
```bash
# lib/infra.sh (extracted from setup-infrastructure.sh)

ensure_static_ip() {
    local ip_name=$1
    local region=$2

    # Check if IP exists
    local existing_ip=$(gcloud compute addresses describe "$ip_name" \
        --region="$region" --format="value(address)" 2>/dev/null || true)

    if [[ -n "$existing_ip" ]]; then
        echo_info "Static IP already exists: $existing_ip"
        echo "$existing_ip"
    else
        echo_info "Creating static IP '$ip_name'..."
        gcloud compute addresses create "$ip_name" --region="$region"
        existing_ip=$(gcloud compute addresses describe "$ip_name" \
            --region="$region" --format="value(address)")
        echo_info "Static IP created: $existing_ip"
        echo "$existing_ip"
    fi
}
```

### Pattern 3: Progressive Feedback

**What:** Even though non-interactive, show progress with timestamped log messages

**When:** Throughout deployment process

**Example:**
```bash
# deploy-auto.sh

echo_info "Starting automated deployment..."
echo_info "[1/6] Loading configuration from config.yaml"
echo_info "[2/6] Finding cheapest region for 64 cores / 256GB RAM"
echo_info "      Checking pricing data for 47 regions..."
echo_info "      Cheapest: us-central1 (n2-standard-64) at $0.78/hour"
echo_info "[3/6] Creating infrastructure (static IP, firewall, disk)"
echo_info "[4/6] Calculating Docker memory limits: 248GB"
echo_info "[5/6] Deploying VM (this takes ~3 minutes)..."
echo_info "[6/6] Waiting for services to start..."
echo_info "Deployment complete! Access at: http://34.44.198.243"
```

### Pattern 4: Atomic Deployment with Rollback

**What:** If deployment fails midway, clean up partial resources

**When:** Any gcloud command fails

**Example:**
```bash
# deploy-auto.sh

set -euo pipefail
trap cleanup ERR

cleanup() {
    echo_error "Deployment failed, cleaning up..."

    # Delete VM if it was created
    if [[ -n "${VM_NAME:-}" ]]; then
        gcloud compute instances delete "$VM_NAME" --zone="$ZONE" --quiet 2>/dev/null || true
    fi

    # Leave infrastructure (static IP, disk) intact for retry
    echo_info "Infrastructure preserved for retry"
    exit 1
}

# Main deployment
load_config
find_cheapest_region
create_infrastructure
deploy_vm  # If this fails, cleanup runs
wait_for_services
```

## Anti-Patterns to Avoid

### Anti-Pattern 1: Interactive Fallbacks

**What:** Adding "fallback to interactive prompts" when config is incomplete

**Why bad:** Defeats the purpose of non-interactive deployment, scripts hang in CI/CD

**Instead:** Fail fast with clear error message about missing config

### Anti-Pattern 2: Hardcoded Pricing

**What:** Embedding price tables directly in script code

**Why bad:** Prices change frequently (spot pricing adjusts daily), quickly becomes inaccurate

**Instead:** Query third-party pricing API, cache for 24 hours, provide flag to skip pricing optimization

### Anti-Pattern 3: Single Monolithic Script

**What:** Putting all logic in one 2000-line `deploy-auto.sh` file

**Why bad:** Untestable, hard to debug, can't reuse components

**Instead:** Modular libraries with single-responsibility functions, orchestrator script sources libraries

### Anti-Pattern 4: Ignoring Existing Infrastructure

**What:** Always creating fresh infrastructure, failing if resources exist

**Why bad:** User loses static IP (breaking DNS), disk data is lost

**Instead:** Idempotent operations - create if missing, reuse if exists

## Scalability Considerations

### At 1 User

**Concerns:** Script reliability, error messages, documentation

**Approach:** Focus on clear validation errors, comprehensive logging

### At 10 Users

**Concerns:** Concurrent deployments (IP address conflicts), shared pricing cache

**Approach:** Unique resource names per deployment, local pricing cache per user

### At 100 Users

**Concerns:** GCP API rate limits, pricing query bottlenecks

**Approach:** Centralized pricing cache service, exponential backoff for gcloud API calls

## Build Order

### Phase 1: Config Parser and Validator
**Goal:** Load YAML config, validate specs, fail fast on bad config

**Deliverables:**
- `lib/config.sh` with `load_config()` and `validate_config()` functions
- `config.yaml.example` with documented schema
- Unit tests for validation logic (use [bats](https://github.com/bats-core/bats-core) or manual test cases)

**Dependencies:** `yq` binary (install from GitHub releases)

### Phase 2: Pricing Query Module
**Goal:** Find cheapest region for given CPU/RAM specs

**Deliverables:**
- `lib/pricing.sh` with `fetch_pricing_data()` and `get_cheapest_region()` functions
- Pricing cache mechanism (store in `/tmp/gcp-pricing-cache-$(date +%Y%m%d).json`)
- Fallback to specific region if pricing query fails

**Dependencies:** CloudPrice.net API or scraping logic, `jq` for JSON parsing

**Known limitation:** Pricing data is MEDIUM confidence, may not reflect real-time spot prices

### Phase 3: Machine Type Mapping
**Goal:** Calculate Docker memory limits and NWCHEM_NPROC from VM specs

**Deliverables:**
- `lib/machine.sh` with `calculate_docker_limits()` function
- Logic to generate `.env` file content
- Startup script template with variable substitution

**Dependencies:** None (pure bash math)

### Phase 4: Orchestrator Script
**Goal:** Non-interactive deployment from config file

**Deliverables:**
- `deploy-auto.sh` that sources all libraries and orchestrates deployment
- Integration with existing infrastructure creation logic
- Progress feedback with timestamps
- Error handling and cleanup on failure

**Dependencies:** All previous phases, existing `lib/infra.sh` and `lib/deploy.sh`

### Phase 5: Testing and Documentation
**Goal:** Validate end-to-end deployment, document usage

**Deliverables:**
- Integration test: deploy with config, verify services start
- README section on automated deployment
- Comparison table: interactive vs. automated deployment
- Troubleshooting guide for common config errors

**Dependencies:** All previous phases, access to GCP project for testing

## Known Issues and Mitigations

### Issue 1: Pricing Data Reliability

**Problem:** Third-party pricing APIs may return stale data or become unavailable

**Mitigation:**
- Cache pricing data with 24-hour TTL
- Provide `--region` flag to override automatic selection
- Fallback to manual pricing query if API fails: `gcloud compute machine-types list --zones | head -n 1`

### Issue 2: Machine Type Availability

**Problem:** Cheapest region may not have spot capacity at deployment time

**Mitigation:**
- Query `gcloud compute zones list` to check for maintenance windows
- Provide top 3 regions ranked by price, try in order until deployment succeeds
- Allow `exclude_regions` config to skip unreliable zones

### Issue 3: Docker nproc Inside Containers

**Problem:** CPU detection inside containers returns 1 (cgroups isolation)

**Mitigation:** Already solved in v2.6 - detect on host, write to `.env` file, use `--oversubscribe` flag

**Source:** `.planning/post-v2.6-problems.md`

### Issue 4: YAML Parser Dependency

**Problem:** `yq` may not be installed on user's system

**Mitigation:**
- Check for `yq` in script, download from GitHub releases if missing
- Alternative: Use simple bash YAML parser (regex-based, LIMITED to flat configs)
- Provide `.env` format as alternative to YAML (simpler, no nested structure)

## Alternative Approaches Considered

### Alternative 1: Use Terraform

**Pros:** Declarative config, state management, mature ecosystem

**Cons:**
- Requires Terraform installation (large dependency)
- Overkill for single-VM deployment
- User must learn Terraform DSL
- Doesn't solve pricing query problem (still need third-party data)

**Verdict:** Not recommended - adds complexity without solving core problems

### Alternative 2: Python Orchestrator

**Pros:** Better data structures, easier API calls, existing libraries (`google-cloud-compute`)

**Cons:**
- Requires Python 3.9+ with pip dependencies
- User already has bash scripts, adding Python is confusing
- Doesn't integrate with existing bash infrastructure scripts

**Verdict:** Not recommended - violates "use existing ecosystem" principle

### Alternative 3: Interactive with Config Override

**Pros:** Keeps existing scripts, adds `--config` flag for automation

**Cons:**
- Half-measures - still need to rip out prompts for CI/CD
- Mixing interactive and non-interactive modes leads to confusing UX
- Doesn't address pricing query requirement

**Verdict:** Not recommended - technical debt accumulation

## Confidence Assessment

| Area | Level | Reason |
|------|-------|--------|
| Config parsing | HIGH | YAML with yq is well-established pattern |
| Pricing query | MEDIUM | Third-party APIs are reliable but not official |
| Machine mapping | HIGH | gcloud machine-types list is official, stable |
| Docker limits | HIGH | Solved in v2.6, pattern is proven |
| Script organization | HIGH | Modular bash is standard practice |

## Summary for Roadmap

**Recommended phase structure:**

1. **Config and validation** - Parse YAML, validate specs, fail fast
2. **Pricing integration** - Query third-party API, find cheapest region
3. **Machine mapping** - Calculate Docker limits, generate startup script
4. **Orchestration** - Tie everything together, refactor existing scripts
5. **Testing and docs** - End-to-end validation, user documentation

**Key architectural decisions:**
- **Modular bash libraries, single entry point** for testability and reuse
- **YAML config with yq parser** for human-friendly configuration
- **Third-party pricing API** (CloudPrice.net) with local caching
- **Host-based CPU/RAM detection** written to `.env` file for Docker
- **Idempotent infrastructure operations** for safe re-runs

**Research flags for phases:**
- Phase 2 (Pricing): MEDIUM confidence on third-party API reliability
- Phase 4 (Orchestration): May need iteration on error handling and rollback logic

## Sources

**GCP Official Documentation:**
- [Get Google Cloud pricing information API](https://docs.cloud.google.com/billing/docs/how-to/get-pricing-information-api)
- [Spot VMs | Compute Engine Documentation](https://docs.cloud.google.com/compute/docs/instances/spot)
- [Pricing | Spot VMs | Google Cloud](https://cloud.google.com/spot-vms/pricing)
- [Use startup scripts on Linux VMs](https://cloud.google.com/compute/docs/instances/startup-scripts/linux)
- [Docker Resource Constraints](https://docs.docker.com/engine/containers/resource_constraints/)

**Third-Party Pricing Tools:**
- [GCP Compute Engine Pricing - CloudPrice](https://cloudprice.net/gcp/compute)
- [GCP Region Pricing Comparison](https://cloudprice.net/gcp/regions)
- [GCP Spot Price History](https://cloudprice.net/gcp/spot-history)

**Docker and Configuration:**
- [Docker Compose Deploy Specification](https://docs.docker.com/reference/compose-file/deploy/)
- [Setting Memory And CPU Limits In Docker | Baeldung](https://www.baeldung.com/ops/docker-memory-limit)
- [Parsing config files with Bash | Opensource.com](https://opensource.com/article/21/6/bash-config)
- [The Ultimate Guide to Modularizing Bash Script Code | Medium](https://medium.com/mkdir-awesome/the-ultimate-guide-to-modularizing-bash-script-code-f4a4d53000c2)

**Best Practices:**
- [How I Built a Production-Ready Deployment Bash Script | Medium](https://medium.com/@ChinecheremU/how-i-built-a-production-ready-deployment-script-for-dockerized-applications-7620f691380d)
- [Bash Script Validation: File Testing, String Checking and Input Validation](https://www.mayhemcode.com/2025/09/bash-script-validation-file-testing.html)
