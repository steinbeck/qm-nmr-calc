# Phase 46: VM Deployment Script - Research

**Researched:** 2026-02-04
**Domain:** GCP Spot VM deployment with Docker Compose
**Confidence:** HIGH

## Summary

This research investigated how to create a single-command GCP Spot VM deployment script that installs Docker, pulls GHCR images, configures HTTPS with Let's Encrypt via Caddy, and handles graceful shutdown during preemption. The phase builds on existing infrastructure scripts from Phase 45 (static IP, firewall rules, persistent disk).

GCP Spot VMs provide 60-91% cost savings over on-demand instances but can be preempted with 30 seconds notice. The standard approach uses `gcloud compute instances create` with `--provisioning-model=SPOT`, startup scripts via `--metadata-from-file` to install Docker and launch containers, and shutdown scripts to gracefully stop containers within the 25-second safe window for stateful workloads. Cost estimation requires querying the Cloud Billing API or displaying machine type specs for user evaluation.

The existing project uses Caddy for automatic HTTPS (configured via DOMAIN environment variable), has Phase 45 infrastructure in place (static IP, firewall, persistent disk), and publishes multi-architecture images to GHCR. The deployment script needs interactive prompts with sensible defaults (region, zone, machine type), cost estimation display, and a Docker Compose override file for GCP-specific settings.

**Primary recommendation:** Use `gcloud compute instances create` with Spot provisioning, attach existing persistent disk, pass startup script via metadata to install Docker Compose and deploy containers, and use shutdown script to trigger `docker compose stop` within 25 seconds.

## Standard Stack

The established tools for GCP Spot VM deployment with Docker:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| gcloud CLI | latest | Create/manage GCP resources | Official GCP command-line tool, industry standard |
| Docker | latest (via apt) | Container runtime | Required for running containerized workloads |
| Docker Compose | v2 (via apt) | Multi-container orchestration | Simplifies deployment of API + worker + Caddy services |
| Caddy | 2.x (alpine image) | Reverse proxy with auto-HTTPS | Automatic Let's Encrypt certificates, zero-config TLS |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Cloud Billing API | v1beta | Programmatic pricing data | Cost estimation feature (optional, can use machine type specs instead) |
| bash | 5.x | Deployment script | Standard for GCP infrastructure automation |
| jq | latest | JSON parsing (optional) | If using Cloud Billing API for cost estimation |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Startup script metadata | cloud-init user-data | cloud-init requires special configuration and is deprecated for container deployment; GCP now recommends startup scripts |
| Spot VMs (--provisioning-model=SPOT) | Preemptible VMs (deprecated flag) | Spot VMs are the current standard; preemptible is legacy |
| Shutdown script | systemd service with preemption detection | systemd approach is more complex; shutdown scripts are GCP's recommended pattern |
| gcloud commands | Terraform | Terraform adds complexity; bash scripts are sufficient for single-VM deployment |

**Installation:**
```bash
# Script installs Docker and Docker Compose via startup script
# On Debian-based images:
apt-get update
apt-get install -y docker.io docker-compose-v2
```

## Architecture Patterns

### Recommended Script Structure
```
gcp/
├── config.sh              # Existing: project, region, zone, prefix, disk size
├── setup-infrastructure.sh # Existing: static IP, firewall, disk (Phase 45)
├── deploy-vm.sh           # NEW: Spot VM creation with prompts
├── startup.sh             # NEW: Installs Docker, starts containers
├── shutdown.sh            # NEW: Graceful container stop on preemption
└── teardown-infrastructure.sh # Existing: cleanup (Phase 45)
```

### Pattern 1: Spot VM Creation with Existing Disk Attachment
**What:** Create Spot VM and attach pre-created persistent disk for data persistence across preemptions
**When to use:** When you need stateful workloads on cost-optimized infrastructure
**Example:**
```bash
# Source: https://docs.cloud.google.com/compute/docs/instances/create-use-spot
# Source: https://docs.cloud.google.com/compute/docs/disks/attach-disks

# Get static IP address
STATIC_IP=$(gcloud compute addresses describe "${RESOURCE_PREFIX}-ip" \
    --region="$GCP_REGION" --format="value(address)")

# Create Spot VM with existing disk attached
gcloud compute instances create "${RESOURCE_PREFIX}-vm" \
    --zone="$GCP_ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --provisioning-model=SPOT \
    --instance-termination-action=STOP \
    --tags="${RESOURCE_PREFIX}-vm" \
    --address="$STATIC_IP" \
    --boot-disk-size=20GB \
    --boot-disk-type=pd-balanced \
    --image-family=debian-12 \
    --image-project=debian-cloud \
    --disk="name=${RESOURCE_PREFIX}-data,mode=rw,device-name=data-disk,boot=no,auto-delete=no" \
    --metadata-from-file=startup-script=startup.sh,shutdown-script=shutdown.sh
```

### Pattern 2: Startup Script for Docker Deployment
**What:** Startup script installs Docker, mounts persistent disk, pulls GHCR images, starts containers
**When to use:** Every Spot VM deployment requiring containerized workloads
**Example:**
```bash
# Source: https://docs.cloud.google.com/compute/docs/instances/startup-scripts/linux
# Source: https://docs.cloud.google.com/compute/docs/containers/deploying-containers

#!/bin/bash
set -euo pipefail

# Install Docker and Docker Compose
apt-get update
apt-get install -y docker.io docker-compose-v2

# Mount persistent disk
DEVICE_NAME="/dev/disk/by-id/google-data-disk"
MOUNT_POINT="/mnt/disks/data"
mkdir -p "$MOUNT_POINT"

# Format if first run (check for existing filesystem)
if ! blkid "$DEVICE_NAME"; then
    mkfs.ext4 -F "$DEVICE_NAME"
fi

# Mount and add to fstab for persistence
mount -o discard,defaults "$DEVICE_NAME" "$MOUNT_POINT"
echo "$DEVICE_NAME $MOUNT_POINT ext4 discard,defaults,nofail 0 2" >> /etc/fstab

# Create project directory and symlink to mounted disk
mkdir -p /opt/qm-nmr-calc
cd /opt/qm-nmr-calc

# Download docker-compose files from project (or embed in metadata)
# Pull images (no authentication needed for public GHCR)
docker compose pull

# Start services
docker compose -f docker-compose.yml -f docker-compose.gcp.yml up -d
```

### Pattern 3: Shutdown Script for Graceful Container Stop
**What:** Shutdown script triggers clean container shutdown within 25-second window during preemption
**When to use:** All Spot VM deployments with stateful containers
**Example:**
```bash
# Source: https://docs.cloud.google.com/compute/docs/shutdownscript
# Source: https://cloud.google.com/kubernetes-engine/docs/concepts/spot-vms

#!/bin/bash
set -euo pipefail

# Change to project directory
cd /opt/qm-nmr-calc || exit 0

# Gracefully stop containers (worker has 300s stop_grace_period internally)
# But preemption only gives us 30s total, aim for 25s to be safe
docker compose stop --timeout 25

# Log completion
echo "$(date): Containers stopped gracefully during preemption" >> /var/log/shutdown.log
```

### Pattern 4: Interactive Prompts with Defaults
**What:** Bash script prompts for region, zone, machine type with sensible defaults
**When to use:** User-facing deployment scripts requiring customization
**Example:**
```bash
# Source: https://ryanstutorials.net/bash-scripting-tutorial/bash-input.php
# Source: https://linuxconfig.org/bash-scripting-how-to-ask-for-user-input

# Load defaults from config.sh
source ./config.sh

# Prompt with defaults
read -p "Select region [${GCP_REGION}]: " INPUT_REGION
SELECTED_REGION="${INPUT_REGION:-${GCP_REGION}}"

read -p "Select zone [${GCP_ZONE}]: " INPUT_ZONE
SELECTED_ZONE="${INPUT_ZONE:-${GCP_ZONE}}"

# List machine types for selected zone
echo "Available machine types in ${SELECTED_ZONE}:"
echo "  e2-standard-2   (2 vCPU, 8 GB RAM)  - Recommended for light workloads"
echo "  e2-standard-4   (4 vCPU, 16 GB RAM) - Recommended for medium workloads"
echo "  n2-standard-4   (4 vCPU, 16 GB RAM) - Higher performance"
echo "  c2-standard-4   (4 vCPU, 16 GB RAM) - Compute-optimized"

read -p "Select machine type [e2-standard-4]: " INPUT_MACHINE
MACHINE_TYPE="${INPUT_MACHINE:-e2-standard-4}"
```

### Pattern 5: Cost Estimation Display
**What:** Display estimated monthly cost before VM creation, allow user to cancel
**When to use:** All cost-sensitive deployments
**Example:**
```bash
# Source: https://docs.cloud.google.com/billing/docs/how-to/get-pricing-information-api
# Simplified approach: show machine specs, let user evaluate cost

echo ""
echo "=========================================="
echo "Cost Estimate for ${MACHINE_TYPE} in ${SELECTED_REGION}"
echo "=========================================="
echo ""
echo "Machine type: ${MACHINE_TYPE}"
echo "Region: ${SELECTED_REGION}"
echo "Provisioning: Spot (60-91% discount vs on-demand)"
echo ""
echo "Estimated monthly cost (Spot pricing, ~730 hours):"
echo "  e2-standard-2:  ~\$15-25/month"
echo "  e2-standard-4:  ~\$30-50/month"
echo "  n2-standard-4:  ~\$50-80/month"
echo "  c2-standard-4:  ~\$70-100/month"
echo ""
echo "Note: Actual costs depend on Spot market pricing and uptime."
echo "Visit https://cloud.google.com/compute/vm-instance-pricing for details."
echo ""
read -p "Continue with VM creation? (yes/no) [yes]: " CONFIRM
CONFIRM="${CONFIRM:-yes}"

if [[ "$CONFIRM" != "yes" ]]; then
    echo "Deployment cancelled."
    exit 0
fi
```

### Pattern 6: Docker Compose GCP Override
**What:** docker-compose.gcp.yml with GCP-specific environment variables and volume paths
**When to use:** When base docker-compose.yml needs cloud environment customization
**Example:**
```yaml
# Source: https://docs.docker.com/compose/how-tos/multiple-compose-files/merge/
# docker-compose.gcp.yml - GCP production overrides

services:
  api:
    volumes:
      # Override to use mounted persistent disk
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info

  caddy:
    environment:
      # Set domain from VM metadata or config
      - DOMAIN=${DOMAIN}

  worker:
    volumes:
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info
      # Let auto-detection handle CPU count on GCP
      - NWCHEM_NPROC=${NWCHEM_NPROC:-}

volumes:
  # Override app-data to not create named volume (use bind mount instead)
  app-data:
    external: false
    driver_opts:
      type: none
      device: /mnt/disks/data
      o: bind
```

### Anti-Patterns to Avoid
- **Attaching disk after VM creation:** Always attach persistent disk during `gcloud compute instances create` using `--disk` flag, not as separate `attach-disk` command. This ensures disk is available for startup script.
- **Using `--address=IP_NAME`:** Must use actual IP address value, not the GCP resource name. Use `gcloud compute addresses describe` to get the IP first.
- **Startup script larger than 256 KB without Cloud Storage:** Scripts over 256 KB must use `startup-script-url` with Cloud Storage. Keep startup scripts small or embed docker-compose files in metadata.
- **Assuming shutdown script always completes:** Shutdown scripts are best-effort. Design containers to handle abrupt termination (e.g., use SIGINT for Huey worker, which is already configured).
- **Using deprecated preemptible flag:** Use `--provisioning-model=SPOT`, not `--preemptible` (deprecated).
- **Querying Cloud Billing API without API key:** Pricing API requires API key setup. For simplicity, display machine type specs and reference pricing URL.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Automatic HTTPS certificates | Custom Let's Encrypt ACME client | Caddy with auto-HTTPS | Caddy handles certificate issuance, renewal, OCSP stapling automatically via `{$DOMAIN}` syntax |
| Container preemption detection | Custom metadata server polling | GCP shutdown scripts with `docker compose stop` | Shutdown scripts are triggered automatically on preemption STOPPING state |
| Disk formatting and mounting | Manual mkfs/mount commands each boot | Startup script with idempotent checks | Use `blkid` to detect existing filesystem, add to `/etc/fstab` for automatic remount |
| Cost estimation API | Cloud Billing API integration | Static cost reference table | API requires authentication setup; hardcoded estimates are simpler and sufficient for deployment script |
| Region/zone selection UI | Custom TUI with dialog/whiptail | Simple `read -p` with defaults | Bash read with defaults is sufficient; fancy UIs add dependency complexity |
| Docker installation | Downloading Docker binaries manually | `apt-get install docker.io docker-compose-v2` | Official Debian packages handle installation, permissions, systemd service setup |

**Key insight:** GCP provides purpose-built primitives (startup/shutdown scripts, metadata, Spot provisioning) that are more reliable than custom solutions. Caddy's automatic HTTPS eliminates the entire certificate management problem.

## Common Pitfalls

### Pitfall 1: Shutdown Script Timeout Exceeding Safe Window
**What goes wrong:** Shutdown script uses `docker compose stop` with container's full `stop_grace_period` (300s for worker), but preemption only allows 30s. Containers are force-killed, risking data corruption.
**Why it happens:** docker-compose inherits stop_grace_period from services, but preemption ignores this. GCP forcibly terminates VM after ~30s.
**How to avoid:** Use `docker compose stop --timeout 25` in shutdown script to override container grace periods and complete within safe window.
**Warning signs:** Logs showing "SIGKILL" instead of "SIGTERM", incomplete job state, corrupted database files after preemption.

### Pitfall 2: Static IP Assignment Using Resource Name Instead of Address
**What goes wrong:** Using `--address="${RESOURCE_PREFIX}-ip"` fails with "Invalid value for field 'address'" error.
**Why it happens:** `--address` flag requires the actual IP address (e.g., "34.123.45.67"), not the GCP resource name.
**How to avoid:** Query the IP address first: `STATIC_IP=$(gcloud compute addresses describe "${RESOURCE_PREFIX}-ip" --region="$GCP_REGION" --format="value(address)")`
**Warning signs:** gcloud error about invalid address format during instance creation.

### Pitfall 3: Persistent Disk Not Mounted Before Docker Compose Starts
**What goes wrong:** Docker creates local volumes instead of using persistent disk data, losing all job history and Let's Encrypt certificates on preemption.
**Why it happens:** Startup script starts `docker compose up` before mounting persistent disk to `/mnt/disks/data`.
**How to avoid:** Mount disk and verify mount point exists before running docker compose. Use bind mount in docker-compose.gcp.yml to override volume paths.
**Warning signs:** Data not persisting across reboots, Let's Encrypt hitting rate limits due to re-issuing certificates.

### Pitfall 4: First Boot Disk Format Errors
**What goes wrong:** `mkfs.ext4` runs on every boot, reformatting disk and destroying data.
**Why it happens:** Startup script doesn't check for existing filesystem before formatting.
**How to avoid:** Use `if ! blkid "$DEVICE_NAME"; then mkfs.ext4 -F "$DEVICE_NAME"; fi` to only format on first boot.
**Warning signs:** Data loss on VM restart, filesystem recreation logs in startup script output.

### Pitfall 5: Caddy Not Getting Valid Let's Encrypt Certificate
**What goes wrong:** Caddy serves self-signed certificate instead of Let's Encrypt, browsers show security warnings.
**Why it happens:** DOMAIN environment variable not set, or DNS not pointing to static IP yet.
**How to avoid:** Verify DNS propagation before deployment (`dig +short $DOMAIN`), ensure DOMAIN set in docker-compose.gcp.yml or .env file.
**Warning signs:** "using self-signed certificate" in Caddy logs, `{$DOMAIN}` resolving to "localhost" in Caddyfile.

### Pitfall 6: Worker Container Calculations Interrupted by Preemption
**What goes wrong:** Long NMR calculations (up to 5 minutes) are lost when Spot VM preempts mid-calculation.
**Why it happens:** Shutdown script stops containers within 25s, but worker has active job running.
**How to avoid:** Accept this limitation for Spot VMs (cost vs completion tradeoff), or use on-demand VMs for critical workloads. Huey worker with SIGINT handles graceful shutdown but can't guarantee completion.
**Warning signs:** Jobs stuck in "started" state after preemption, need manual retry.

### Pitfall 7: Docker Not Available in PATH for Startup Script
**What goes wrong:** Startup script fails at `docker compose pull` with "command not found".
**Why it happens:** Docker installed via apt but systemd service not started yet when startup script runs.
**How to avoid:** Add `systemctl enable docker && systemctl start docker` after installation, or use `sleep 5` to wait for service startup.
**Warning signs:** Startup script errors about docker command not found, containers not running after VM creation.

## Code Examples

Verified patterns from official sources:

### Creating Spot VM with All Required Flags
```bash
# Source: https://docs.cloud.google.com/compute/docs/instances/create-use-spot
# Source: https://docs.cloud.google.com/compute/docs/disks/attach-disks

#!/bin/bash
set -euo pipefail

# Get static IP address value (not name)
STATIC_IP=$(gcloud compute addresses describe "${RESOURCE_PREFIX}-ip" \
    --region="$GCP_REGION" \
    --format="value(address)")

# Create Spot VM
gcloud compute instances create "${RESOURCE_PREFIX}-vm" \
    --zone="$GCP_ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --provisioning-model=SPOT \
    --instance-termination-action=STOP \
    --network-interface=address="$STATIC_IP" \
    --tags="${RESOURCE_PREFIX}-vm" \
    --boot-disk-size=20GB \
    --boot-disk-type=pd-balanced \
    --image-family=debian-12 \
    --image-project=debian-cloud \
    --disk="name=${RESOURCE_PREFIX}-data,mode=rw,device-name=data-disk,boot=no,auto-delete=no" \
    --metadata-from-file=startup-script=startup.sh,shutdown-script=shutdown.sh \
    --metadata=DOMAIN="$DOMAIN"
```

### Complete Startup Script with Docker Installation
```bash
# Source: https://docs.cloud.google.com/compute/docs/instances/startup-scripts/linux
# Source: https://docs.cloud.google.com/compute/docs/containers/deploying-containers

#!/bin/bash
set -euo pipefail

# Log startup script execution
exec > >(tee -a /var/log/startup-script.log)
exec 2>&1
echo "$(date): Starting VM setup..."

# Install Docker and Docker Compose
echo "Installing Docker..."
apt-get update
apt-get install -y docker.io docker-compose-v2 curl

# Ensure Docker service is running
systemctl enable docker
systemctl start docker

# Wait for Docker to be ready
for i in {1..10}; do
    if docker info >/dev/null 2>&1; then
        echo "Docker is ready"
        break
    fi
    echo "Waiting for Docker... ($i/10)"
    sleep 2
done

# Mount persistent disk
DEVICE_NAME="/dev/disk/by-id/google-data-disk"
MOUNT_POINT="/mnt/disks/data"
echo "Mounting persistent disk..."

mkdir -p "$MOUNT_POINT"

# Format only if first boot
if ! blkid "$DEVICE_NAME"; then
    echo "Formatting disk (first boot)..."
    mkfs.ext4 -F "$DEVICE_NAME"
fi

# Mount disk
mount -o discard,defaults "$DEVICE_NAME" "$MOUNT_POINT"

# Add to fstab for automatic remount on restart
if ! grep -q "$DEVICE_NAME" /etc/fstab; then
    echo "$DEVICE_NAME $MOUNT_POINT ext4 discard,defaults,nofail 0 2" >> /etc/fstab
fi

# Create application directory
mkdir -p /opt/qm-nmr-calc
cd /opt/qm-nmr-calc

# Fetch DOMAIN from instance metadata
DOMAIN=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/attributes/DOMAIN" -H "Metadata-Flavor: Google")
echo "DOMAIN=${DOMAIN}" > .env

# Create docker-compose files (inline or download from Cloud Storage)
# For this example, assume files are embedded in metadata or downloaded
# In production, use --metadata-from-file to embed or download from GCS

cat > docker-compose.gcp.yml <<'EOF'
services:
  api:
    volumes:
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production

  caddy:
    environment:
      - DOMAIN=${DOMAIN}

  worker:
    volumes:
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production
EOF

# Pull images from GHCR (public, no auth needed)
echo "Pulling images from GHCR..."
docker compose pull

# Start services
echo "Starting services..."
docker compose -f docker-compose.yml -f docker-compose.gcp.yml up -d

echo "$(date): VM setup complete!"
```

### Complete Shutdown Script for Graceful Preemption
```bash
# Source: https://docs.cloud.google.com/compute/docs/shutdownscript
# Source: https://cloud.google.com/kubernetes-engine/docs/concepts/spot-vms

#!/bin/bash
set -euo pipefail

# Log shutdown script execution
exec > >(tee -a /var/log/shutdown-script.log)
exec 2>&1
echo "$(date): Preemption detected, stopping containers..."

# Change to project directory
cd /opt/qm-nmr-calc || exit 0

# Gracefully stop containers within safe window (25s target)
# This overrides container stop_grace_period settings
if docker compose ps -q 2>/dev/null; then
    echo "Stopping Docker Compose services..."
    docker compose stop --timeout 25
    echo "Containers stopped gracefully"
else
    echo "No containers running, skipping stop"
fi

echo "$(date): Shutdown script complete"
```

### Interactive Prompt Pattern with Validation
```bash
# Source: https://linuxconfig.org/bash-scripting-how-to-ask-for-user-input

#!/bin/bash
set -euo pipefail

# Load defaults from config
source ./config.sh

# Region selection with validation
while true; do
    read -p "Select region [${GCP_REGION}]: " INPUT_REGION
    SELECTED_REGION="${INPUT_REGION:-${GCP_REGION}}"

    # Validate region exists
    if gcloud compute regions describe "$SELECTED_REGION" &>/dev/null; then
        break
    else
        echo "Error: Region '$SELECTED_REGION' not found. Try 'us-central1', 'europe-west1', etc."
    fi
done

# Zone selection with validation
while true; do
    read -p "Select zone [${GCP_ZONE}]: " INPUT_ZONE
    SELECTED_ZONE="${INPUT_ZONE:-${GCP_ZONE}}"

    if gcloud compute zones describe "$SELECTED_ZONE" &>/dev/null; then
        break
    else
        echo "Error: Zone '$SELECTED_ZONE' not found. Try '${SELECTED_REGION}-a', '${SELECTED_REGION}-b', etc."
    fi
done

# Machine type with helper information
echo ""
echo "Available machine types (Spot discount ~60-91%):"
echo "  e2-standard-2   (2 vCPU, 8 GB)   ~\$15-25/month  - Light workloads"
echo "  e2-standard-4   (4 vCPU, 16 GB)  ~\$30-50/month  - Recommended"
echo "  n2-standard-4   (4 vCPU, 16 GB)  ~\$50-80/month  - Higher performance"
echo "  c2-standard-4   (4 vCPU, 16 GB)  ~\$70-100/month - Compute-optimized"
echo ""
read -p "Select machine type [e2-standard-4]: " INPUT_MACHINE
MACHINE_TYPE="${INPUT_MACHINE:-e2-standard-4}"

# Validate machine type exists in zone
if ! gcloud compute machine-types describe "$MACHINE_TYPE" --zone="$SELECTED_ZONE" &>/dev/null; then
    echo "Warning: Machine type '$MACHINE_TYPE' not available in zone '$SELECTED_ZONE'"
    echo "Defaulting to e2-standard-4"
    MACHINE_TYPE="e2-standard-4"
fi
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Preemptible VMs (--preemptible flag) | Spot VMs (--provisioning-model=SPOT) | 2021-2022 | Spot VMs are the new standard; preemptible flag still works but deprecated |
| GCP built-in container deployment (konlet) | Startup scripts with Docker commands | Deprecated 2023 | GCP now recommends startup scripts; konlet feature removed from documentation |
| cloud-init for container startup | Metadata startup scripts | Ongoing preference | Startup scripts are simpler and officially recommended over cloud-init |
| Manual Let's Encrypt with certbot | Caddy automatic HTTPS | Caddy 2.x (2020+) | Eliminates certificate management complexity entirely |
| Docker Compose v1 (docker-compose) | Docker Compose v2 (docker compose) | 2022-2023 | V2 is standalone binary, better performance, integrated with Docker CLI |
| Separate docker instances attach-disk command | --disk flag in instances create | Always available | Single command is more reliable; disk ready for startup script |

**Deprecated/outdated:**
- **konlet (GCP container deployment feature)**: Deprecated; GCP documentation removed this feature. Use startup scripts with `docker run` or `docker compose` instead.
- **--preemptible flag**: Still works but deprecated. Use `--provisioning-model=SPOT` instead.
- **startup-script-url with gs://**: Still works but startup-script (inline or local file) is simpler for scripts under 256 KB.
- **Docker Compose v1 syntax**: docker-compose (hyphenated) is legacy. Use docker compose (space) with v2 plugin.

## Open Questions

Things that couldn't be fully resolved:

1. **Cloud Billing API practical usability**
   - What we know: API exists for querying VM pricing programmatically; requires API key setup
   - What's unclear: Whether API returns Spot pricing (dynamic) or just on-demand pricing; API response structure for machine types
   - Recommendation: Start with hardcoded cost estimates in script; add API integration as enhancement if users request it

2. **Docker Compose file distribution strategy**
   - What we know: Startup scripts can embed small files or download from Cloud Storage
   - What's unclear: Best practice for distributing docker-compose.yml + docker-compose.gcp.yml to VM (embed in metadata, download from GCS bucket, git clone repo)
   - Recommendation: Embed docker-compose.gcp.yml in startup script; reference docker-compose.yml from GHCR release or git repo. Keep startup script under 256 KB to avoid Cloud Storage requirement.

3. **Optimal machine type default for NMR calculations**
   - What we know: Worker uses NWChem with MPI (CPU-bound), auto-detects CPU count
   - What's unclear: Whether e2-standard-4 (4 vCPU, 16 GB) is sufficient or if c2-standard-4 (compute-optimized) provides better performance/cost ratio
   - Recommendation: Default to e2-standard-4 for cost-effectiveness; document c2-standard-4 as high-performance option. User can test and reconfigure.

4. **DNS propagation timing during deployment**
   - What we know: Caddy needs DNS pointing to static IP before Let's Encrypt certificate issuance succeeds
   - What's unclear: Whether script should check DNS propagation or just document prerequisite
   - Recommendation: Document DNS as prerequisite in script output; optionally add `dig +short $DOMAIN` check with warning if not propagated yet

5. **Handling worker jobs during preemption**
   - What we know: Worker has 5-minute stop_grace_period, but preemption forces 25s shutdown
   - What's unclear: Whether Huey worker checkpoints job state or loses progress entirely
   - Recommendation: Document that in-progress jobs may be lost during preemption; Huey should retry failed jobs. This is acceptable tradeoff for 60-91% cost savings.

## Sources

### Primary (HIGH confidence)
- [GCP Spot VMs Documentation](https://docs.cloud.google.com/compute/docs/instances/spot) - Preemption process, 30-second window, shutdown scripts
- [Create and Use Spot VMs](https://docs.cloud.google.com/compute/docs/instances/create-use-spot) - gcloud command syntax, --provisioning-model=SPOT flag
- [Startup Scripts on Linux VMs](https://docs.cloud.google.com/compute/docs/instances/startup-scripts/linux) - Metadata keys, execution timing, best practices
- [Run Shutdown Scripts](https://docs.cloud.google.com/compute/docs/shutdownscript) - Shutdown script metadata, execution timing, limitations
- [Attach Disks to VMs](https://docs.cloud.google.com/compute/docs/disks/attach-disks) - --disk flag syntax for attaching existing disks
- [GKE Spot VMs](https://cloud.google.com/kubernetes-engine/docs/concepts/spot-vms) - 25-second graceful termination recommendation for stateful workloads
- [Cloud Billing API](https://docs.cloud.google.com/billing/docs/how-to/get-pricing-information-api) - Pricing API endpoints and authentication
- [Caddy Automatic HTTPS](https://caddyserver.com/docs/automatic-https) - Let's Encrypt integration, ACME protocol
- [Docker Compose Merge](https://docs.docker.com/compose/how-tos/multiple-compose-files/merge/) - Override file pattern for environment-specific settings
- [Configure Static External IP Addresses](https://docs.cloud.google.com/compute/docs/ip-addresses/configure-static-external-ip-address) - Static IP assignment with --address flag
- [Add Network Tags](https://docs.cloud.google.com/vpc/docs/add-remove-network-tags) - --tags flag for firewall rule targeting

### Secondary (MEDIUM confidence)
- [GCP Spot VMs Explained: A Smarter Way to Cut Cloud Costs](https://www.pump.co/blog/spot-instances-gcp) - 60-91% discount figures, practical usage patterns
- [Deploying Containers on Instances](https://cloud.google.com/compute/docs/containers/deploying-containers) - Note that built-in container feature is deprecated
- [How to Create Docker Compose Override Strategies](https://oneuptime.com/blog/post/2026-01-30-docker-compose-override-strategies/view) - Override file best practices (Jan 2026)
- [Bash Scripting: Read User Input](https://linuxconfig.org/bash-scripting-how-to-ask-for-user-input) - Interactive prompt patterns with defaults

### Tertiary (LOW confidence)
- [Run docker-compose on GCE](https://jedri.medium.com/run-docker-compose-on-gce-a5a49f64225e) - Community example of Docker Compose on GCE
- [Reliably Executing Shutdown Scripts in GCE](https://haggainuchi.com/shutdown.html) - Community insights on shutdown script reliability

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official GCP documentation and established tools (gcloud, Docker, Caddy)
- Architecture: HIGH - Verified patterns from official GCP documentation, existing Phase 45 scripts provide foundation
- Pitfalls: HIGH - Pitfalls derived from official documentation limitations and common deployment errors
- Cost estimation: MEDIUM - Pricing figures from community sources; Cloud Billing API structure unclear without testing

**Research date:** 2026-02-04
**Valid until:** 2026-03-04 (30 days - GCP stable platform, pricing updates quarterly)
