# Phase 52: HTTP-Only Container Deployment - Research

**Researched:** 2026-02-06
**Domain:** Docker Compose configuration management, GCP VM lifecycle scripts, HTTP-only deployment
**Confidence:** HIGH

## Summary

Phase 52 focuses on verifying HTTP-only container deployment and ensuring lifecycle scripts work with v2.7's automated deployment system. The critical insight is that Phase 50's `generate_startup_script()` **already creates** the docker-compose.gcp.yml override and .env file, so this phase is primarily about validation and script migration rather than implementing container deployment from scratch.

The main work is migrating v2.6's interactive lifecycle scripts (which use config.sh with hardcoded zone/region) to work with v2.7's TOML-based config system where zone/region are dynamically selected at deployment time. This requires either storing deployment metadata or querying GCP to discover the VM's zone at runtime.

**Key findings:**
- Startup script generation is complete (Phase 50) - HTTP-only docker-compose.gcp.yml override already implemented
- Lifecycle scripts currently source gcp/config.sh (v2.6 bash variables) which is incompatible with v2.7's TOML config + dynamic zone selection
- GCP metadata server provides runtime zone detection without config files: `curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone`
- Teardown script already removes all infrastructure (IP, firewall, disk) - just needs HTTPS rule removal verification

**Primary recommendation:** Use GCP metadata server for runtime zone detection in lifecycle scripts, eliminating dependency on static config files while preserving v2.6 script functionality.

## Standard Stack

The established tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Docker Compose | 2.x | Container orchestration | Industry standard for multi-container applications |
| gcloud CLI | Latest | GCP compute management | Official Google Cloud SDK, required for all GCP operations |
| GCP Metadata Server | N/A | Runtime VM information | Built-in to all GCP VMs, authoritative source for zone/region |
| bash | 4.x+ | Shell scripting | Universal Unix shell, installed by default on GCP VMs |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| jq | 1.6+ | JSON parsing in bash | Extract fields from gcloud JSON output |
| curl | Latest | HTTP requests | Query metadata server, health checks |
| Python 3.11+ | 3.11+ | Config validation | Already used for validate_config.py and select_machine.py |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| GCP Metadata Server | Store zone in file | Metadata server is authoritative and always available |
| bash scripts | Python CLI | Bash scripts are simpler for gcloud wrappers, no new dependencies |
| TOML config for lifecycle | Extend v2.6 config.sh | TOML config doesn't help when zone is dynamic |

**Installation:**
All tools are pre-installed on GCP VMs (gcloud, curl, bash). jq is installed by startup script during deployment.

## Architecture Patterns

### Recommended Project Structure
```
gcp/
├── lib/
│   ├── config.sh           # TOML config loader (Phase 49)
│   ├── pricing.sh          # Pricing query library (Phase 49)
│   ├── machine.sh          # Machine selection wrapper (Phase 50)
│   └── infra.sh            # Infrastructure operations (Phase 51)
├── deploy-auto.sh          # Main deployment orchestrator (Phase 51)
├── start-vm.sh             # Start stopped VM (migrate for v2.7)
├── stop-vm.sh              # Stop running VM (migrate for v2.7)
├── delete-vm.sh            # Delete VM, preserve disk (migrate for v2.7)
├── status-vm.sh            # Show VM status (migrate for v2.7)
├── ssh-vm.sh               # SSH into VM (migrate for v2.7)
├── logs-vm.sh              # View container logs (migrate for v2.7)
├── teardown-infrastructure.sh  # Remove all resources (verify HTTPS cleanup)
├── select_machine.py       # Machine selection module (Phase 50)
├── validate_config.py      # Config validation (Phase 49)
└── query_pricing.py        # Pricing query (Phase 49)
```

### Pattern 1: GCP Metadata Server Zone Detection
**What:** Query VM's zone at runtime from GCP metadata server instead of static config files
**When to use:** In lifecycle scripts (start, stop, delete, status, ssh, logs) that need to know which zone the VM is in
**Example:**
```bash
# Source: https://cloud.google.com/compute/docs/metadata/querying-metadata
get_vm_zone() {
    local vm_name="${1}"

    # Try metadata server first (if running on GCP VM)
    local zone
    zone=$(curl -sf -H "Metadata-Flavor: Google" \
        http://metadata.google.internal/computeMetadata/v1/instance/zone 2>/dev/null | cut -d/ -f4)

    if [[ -n "$zone" ]]; then
        echo "$zone"
        return 0
    fi

    # Fallback: Query gcloud for VM location (works from anywhere)
    gcloud compute instances list \
        --filter="name=${vm_name}" \
        --format="value(zone)" 2>/dev/null | head -1
}
```

### Pattern 2: HTTP-Only Docker Compose Override
**What:** Override base docker-compose.yml to expose HTTP on port 80 without Caddy/HTTPS services
**When to use:** GCP spot instances where HTTPS/domain configuration is not desired (v2.7 requirement)
**Example:**
```yaml
# Source: Generated by gcp/select_machine.py generate_startup_script()
# docker-compose.gcp.yml
services:
  api:
    volumes:
      - /mnt/disks/data:/app/data  # Use mounted persistent disk
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info
    ports:
      - "80:8000"  # HTTP only, map to host port 80

  worker:
    volumes:
      - /mnt/disks/data:/app/data  # Use mounted persistent disk
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info
    deploy:
      resources:
        limits:
          memory: ${WORKER_MEMORY_LIMIT}  # Dynamically calculated

# Caddy service completely omitted (HTTP-only deployment)
volumes:
  app-data: null  # Disable named volume, use bind mount instead
```

### Pattern 3: Container Health Check Validation
**What:** Verify containers are healthy after startup using Docker Compose health checks
**When to use:** After VM creation or restart to confirm deployment succeeded
**Example:**
```bash
# Source: https://docs.docker.com/compose/how-tos/startup-order/
wait_for_healthy() {
    local max_wait=180  # 3 minutes
    local elapsed=0

    echo "Waiting for containers to become healthy..."

    while [[ $elapsed -lt $max_wait ]]; do
        # Check health status of all services with health checks
        local unhealthy
        unhealthy=$(docker compose ps --format json | \
            jq -r '.[] | select(.Health == "unhealthy" or .Health == "starting") | .Service' | \
            wc -l)

        if [[ "$unhealthy" -eq 0 ]]; then
            echo "All containers healthy!"
            return 0
        fi

        sleep 5
        elapsed=$((elapsed + 5))
    done

    echo "ERROR: Containers did not become healthy within ${max_wait}s" >&2
    docker compose ps
    return 1
}
```

### Anti-Patterns to Avoid
- **Hardcoding zone/region in lifecycle scripts:** v2.7 selects zone dynamically at deployment time, hardcoded values will be wrong for automated deployments
- **Storing deployment metadata in files:** Files can be lost or desync'd, metadata server is always authoritative
- **Checking container status with docker ps:** Use `docker compose ps` with health checks for accurate service status
- **Assuming containers are ready when started:** Docker Compose starts containers quickly but apps need time to initialize - use health checks

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| VM zone detection | Parse gcloud output manually | GCP Metadata Server API | Authoritative, no parsing, works inside VM |
| Container readiness | Poll with curl loops | Docker Compose health checks | Built-in, standardized, handles retries/timeouts |
| JSON parsing in bash | awk/sed/grep | jq | Handles nested objects, arrays, escaping correctly |
| Dynamic Docker Compose overrides | Script-generated YAML | Environment variables + override file | Compose merging is battle-tested, handles edge cases |
| Script error handling | Check $? manually | set -euo pipefail | Catches errors in pipes, undefined variables, subshells |

**Key insight:** GCP provides built-in infrastructure for VM metadata queries and Docker Compose has mature health check patterns. Using these instead of custom solutions reduces bugs and maintenance burden.

## Common Pitfalls

### Pitfall 1: Assuming Static Zone Configuration
**What goes wrong:** Lifecycle scripts source config.sh expecting GCP_ZONE to be set, but v2.7 selects zone dynamically at deployment time based on pricing
**Why it happens:** v2.6 used interactive prompts with hardcoded zone in config.sh, v2.7 automated this away
**How to avoid:** Query GCP for the VM's actual zone at runtime using either gcloud or metadata server
**Warning signs:** Scripts fail with "VM not found in zone us-central1-a" when VM was actually created in europe-west1-b

### Pitfall 2: HTTPS Firewall Rules in Teardown
**What goes wrong:** Teardown script tries to delete ${RESOURCE_PREFIX}-allow-https firewall rule, but v2.7 only creates HTTP and SSH rules (no HTTPS)
**Why it happens:** v2.6 created HTTPS rules for Caddy, v2.7 is HTTP-only
**How to avoid:** Make firewall rule deletion idempotent - check if rule exists before trying to delete, or ignore deletion failures
**Warning signs:** Teardown script shows warnings about missing HTTPS rule every time it runs

### Pitfall 3: Hardcoded Caddy References
**What goes wrong:** logs-vm.sh offers "caddy" as a service option, but docker-compose.gcp.yml doesn't include Caddy
**Why it happens:** Script was written for v2.6 HTTPS deployment which included Caddy reverse proxy
**How to avoid:** Remove Caddy-specific options from lifecycle scripts for v2.7
**Warning signs:** `./logs-vm.sh caddy` fails with "no such service: caddy"

### Pitfall 4: Port 80 Binding Conflicts
**What goes wrong:** API container can't bind to port 80 if another service (like Apache) is already using it
**Why it happens:** Some GCP VM images come with pre-installed web servers
**How to avoid:** Startup script should check for port 80 usage and stop conflicting services before starting containers
**Warning signs:** "bind: address already in use" in container logs

### Pitfall 5: Container Memory Limits Too Low
**What goes wrong:** Worker container gets OOM-killed during large NMR calculations because WORKER_MEMORY_LIMIT is too low
**Why it happens:** Docker limit formula (VM_RAM - 8GB) leaves insufficient margin for calculation peaks
**How to avoid:** Phase 50's calculate_docker_resources() validates minimum 4GB after overhead, but monitor actual usage
**Warning signs:** Worker container restarts mid-calculation with exit code 137 (OOM kill)

## Code Examples

Verified patterns from implementation:

### Zone Detection with Metadata Server
```bash
# Source: https://cloud.google.com/compute/docs/metadata/querying-metadata
get_vm_zone() {
    local vm_name="${1}"
    local project_id="${2:-$(gcloud config get-value project 2>/dev/null)}"

    # Query gcloud for VM's zone (works from anywhere with gcloud auth)
    gcloud compute instances list \
        --project="$project_id" \
        --filter="name=${vm_name}" \
        --format="value(zone)" \
        --quiet 2>/dev/null | head -1
}

# Usage in lifecycle scripts
VM_NAME="${RESOURCE_PREFIX}-vm"
GCP_ZONE=$(get_vm_zone "$VM_NAME") || {
    echo "ERROR: Could not determine zone for VM $VM_NAME" >&2
    exit 1
}
```

### HTTP-Only Container Startup Validation
```bash
# Source: Generated by Phase 50 startup script
validate_deployment() {
    local max_wait=180
    local elapsed=0

    echo "Waiting for containers to start..."

    while [[ $elapsed -lt $max_wait ]]; do
        # Check if docker compose services are running
        if docker compose -f /opt/qm-nmr-calc/docker-compose.yml \
           -f /opt/qm-nmr-calc/docker-compose.gcp.yml ps --quiet api worker | grep -q .; then
            echo "Containers started, checking health..."

            # Wait for health checks to pass
            sleep 10
            local unhealthy
            unhealthy=$(docker compose -f /opt/qm-nmr-calc/docker-compose.yml \
                -f /opt/qm-nmr-calc/docker-compose.gcp.yml ps --format json | \
                jq -r '.[] | select(.Health == "unhealthy") | .Service')

            if [[ -z "$unhealthy" ]]; then
                echo "✓ Deployment healthy"
                return 0
            fi
        fi

        sleep 5
        elapsed=$((elapsed + 5))
    done

    echo "ERROR: Deployment validation timed out" >&2
    docker compose ps
    return 1
}
```

### Idempotent Firewall Rule Deletion
```bash
# Source: Pattern from existing teardown-infrastructure.sh
delete_firewall_rule_safe() {
    local rule_name="${1}"

    if gcloud compute firewall-rules describe "$rule_name" &>/dev/null; then
        echo "Deleting firewall rule '$rule_name'..."
        gcloud compute firewall-rules delete "$rule_name" --quiet || {
            echo "WARNING: Failed to delete firewall rule '$rule_name' (may not exist)" >&2
        }
    else
        echo "Firewall rule '$rule_name' does not exist, skipping"
    fi
}

# Usage in teardown script
HTTP_RULE="${RESOURCE_PREFIX}-allow-http"
SSH_RULE="${RESOURCE_PREFIX}-allow-ssh"
# Note: No HTTPS rule for v2.7 HTTP-only deployment

delete_firewall_rule_safe "$HTTP_RULE"
delete_firewall_rule_safe "$SSH_RULE"
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| config.sh with hardcoded zone | Dynamic zone selection + runtime detection | v2.7 (Phase 50) | Lifecycle scripts must query GCP for zone |
| Caddy + HTTPS on GCP | HTTP-only on port 80 | v2.7 (Phase 52) | No Caddy service in docker-compose.gcp.yml |
| Interactive deploy-vm.sh | Non-interactive deploy-auto.sh | v2.7 (Phase 51) | Lifecycle scripts remain interactive |
| Manual docker-compose.gcp.yml | Generated by startup script | v2.7 (Phase 50) | Guaranteed correct resource limits |
| Static NWCHEM_NPROC | Runtime $(nproc) detection | v2.6 fix | Adapts to actual VM CPU count |

**Deprecated/outdated:**
- **gcp/config.sh (v2.6 bash config):** Replaced by gcp/lib/config.sh (TOML loader) for deployment, but lifecycle scripts need different approach (runtime zone detection)
- **HTTPS firewall rules:** v2.7 creates HTTP and SSH only, no HTTPS rule to delete
- **Caddy service references:** docker-compose.gcp.yml omits Caddy entirely for HTTP-only deployment
- **Domain metadata on VM:** v2.7 doesn't use domains, don't check for DOMAIN metadata key

## Open Questions

Things that couldn't be fully resolved:

1. **Should lifecycle scripts support both v2.6 and v2.7 deployments?**
   - What we know: v2.6 VMs have zone in config.sh, v2.7 VMs don't
   - What's unclear: Whether users might have both v2.6 and v2.7 VMs running simultaneously
   - Recommendation: Detect v2.7 VMs by absence of config.sh, fall back to gcloud zone query (works for both)

2. **How to handle zone detection failures?**
   - What we know: gcloud query works from anywhere with auth, metadata server only works inside VM
   - What's unclear: Should scripts try multiple detection methods or fail fast?
   - Recommendation: Single method (gcloud query) with clear error message is simpler and works everywhere

3. **Should status-vm.sh show HTTP URL?**
   - What we know: v2.7 is HTTP-only, status script currently shows "Web UI: https://$DOMAIN" for v2.6
   - What's unclear: Should it show "http://$EXTERNAL_IP" for v2.7? Users might not expect HTTP URL
   - Recommendation: Show "HTTP endpoint: http://$EXTERNAL_IP" to make it clear this is HTTP, not HTTPS

## Sources

### Primary (HIGH confidence)
- GCP Compute Engine documentation (https://cloud.google.com/compute/docs/metadata/querying-metadata) - VM metadata server API
- GCP gcloud SDK documentation (https://cloud.google.com/compute/docs/instances/stop-start-instance) - Instance lifecycle management
- Docker Compose documentation (https://docs.docker.com/compose/how-tos/startup-order/) - Health checks and startup ordering
- Project codebase: gcp/select_machine.py, gcp/lib/config.sh, deploy-auto.sh, existing lifecycle scripts

### Secondary (MEDIUM confidence)
- [Docker Compose port override patterns](https://www.codestudy.net/blog/docker-compose-override-a-ports-property-instead-of-merging-it/) - Port mapping replacement in override files
- [GCP startup and shutdown scripts guide](https://medium.com/google-cloud/few-tips-and-tricks-with-gce-startup-script-323433e2b5ee) - Startup script best practices

### Tertiary (LOW confidence)
- None - all findings verified with official documentation or project codebase

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All tools are standard GCP/Docker tooling with official documentation
- Architecture: HIGH - Patterns verified in existing Phase 50/51 implementation and official GCP docs
- Pitfalls: HIGH - Based on v2.6 → v2.7 migration issues documented in STATE.md and code inspection

**Research date:** 2026-02-06
**Valid until:** 2026-03-06 (30 days - stable technologies, official GCP APIs unlikely to change)
