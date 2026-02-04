# Architecture: GCP Spot VM Deployment

**Project:** qm-nmr-calc GCP Spot Deployment
**Researched:** 2026-02-04
**Focus:** Integrating GCP Spot VMs with existing Docker Compose architecture

## Executive Summary

This research defines how GCP Spot VM deployment integrates with the existing qm-nmr-calc Docker Compose architecture. The key insight is that the existing architecture is already well-suited for spot deployment with minimal modifications. The primary additions are: (1) GCP-specific startup/shutdown scripts, (2) persistent disk for data survival across preemption, and (3) graceful shutdown handling for the worker container.

**Core Principle:** The Docker Compose stack remains unchanged. GCP infrastructure wraps around it.

---

## Integration Architecture

### High-Level Overview

```
+------------------------------------------------------------------+
|                        GCP Spot VM                                |
|                                                                   |
|  +------------------------------------------------------------+  |
|  |                   Startup Script                            |  |
|  |  1. Mount persistent disk to /mnt/data                     |  |
|  |  2. Install Docker (if needed)                             |  |
|  |  3. Clone/pull repo                                        |  |
|  |  4. docker compose up -d                                   |  |
|  +------------------------------------------------------------+  |
|                              |                                    |
|  +------------------------------------------------------------+  |
|  |                  Docker Compose Stack                       |  |
|  |  (UNCHANGED from existing architecture)                     |  |
|  |                                                            |  |
|  |  +--------+  +--------+  +--------+                        |  |
|  |  | Caddy  |  |  API   |  | Worker |                        |  |
|  |  | :80/443|  | :8000  |  | NWChem |                        |  |
|  |  +--------+  +--------+  +--------+                        |  |
|  |       |           |           |                             |  |
|  |       +-----+-----+-----+-----+                             |  |
|  |             |           |                                   |  |
|  |      [app-data]   [caddy_data]                             |  |
|  |             |           |                                   |  |
|  +-------------|-----------|----------------------------------+  |
|                |           |                                     |
|  +-------------+-----------+----------------------------------+  |
|  |              Persistent Disk (/mnt/data)                   |  |
|  |  - Docker volumes bind-mounted here                        |  |
|  |  - Survives preemption/restart                            |  |
|  +------------------------------------------------------------+  |
|                                                                   |
|  +------------------------------------------------------------+  |
|  |                   Shutdown Script                           |  |
|  |  1. docker compose down (graceful, up to 30s)             |  |
|  |  2. Sync filesystems                                       |  |
|  +------------------------------------------------------------+  |
+------------------------------------------------------------------+
```

### Component Interaction

| Component | Location | Purpose | Modified? |
|-----------|----------|---------|-----------|
| Docker Compose | VM | Container orchestration | NO - existing file |
| Caddy | Container | HTTPS termination | NO |
| API | Container | FastAPI service | NO |
| Worker | Container | NWChem calculations | NO |
| Startup script | GCP metadata | Bootstrap on boot | NEW |
| Shutdown script | GCP metadata | Graceful stop | NEW |
| Persistent disk | GCP | Data persistence | NEW |

---

## GCP-Specific Components

### 1. Startup Script

The startup script runs on every VM boot, including after preemption recovery.

**Metadata key:** `startup-script`

**Script responsibilities:**
1. Mount the persistent disk (if not already in fstab)
2. Ensure Docker is installed
3. Pull/update the application
4. Start Docker Compose services

**Example startup script:**

```bash
#!/bin/bash
set -e

# Configuration
DISK_DEVICE="/dev/sdb"
MOUNT_POINT="/mnt/data"
APP_DIR="/opt/qm-nmr-calc"
REPO_URL="https://github.com/steinbeck/qm-nmr-calc.git"

# Log to serial console for debugging
exec > >(tee /var/log/startup-script.log) 2>&1
echo "=== Startup script started at $(date) ==="

# 1. Mount persistent disk if not already mounted
if ! mountpoint -q "$MOUNT_POINT"; then
    echo "Mounting persistent disk..."
    mkdir -p "$MOUNT_POINT"

    # Check if disk needs formatting (first boot only)
    if ! blkid "$DISK_DEVICE" | grep -q ext4; then
        echo "Formatting disk with ext4..."
        mkfs.ext4 -F "$DISK_DEVICE"
    fi

    mount "$DISK_DEVICE" "$MOUNT_POINT"

    # Add to fstab for future boots
    if ! grep -q "$DISK_DEVICE" /etc/fstab; then
        echo "$DISK_DEVICE $MOUNT_POINT ext4 defaults,nofail 0 2" >> /etc/fstab
    fi
fi

# Create data directories on persistent disk
mkdir -p "$MOUNT_POINT/docker-volumes/app-data"
mkdir -p "$MOUNT_POINT/docker-volumes/caddy-data"
mkdir -p "$MOUNT_POINT/docker-volumes/caddy-config"

# 2. Install Docker if not present
if ! command -v docker &> /dev/null; then
    echo "Installing Docker..."
    curl -fsSL https://get.docker.com | sh
fi

# 3. Clone or update repo
if [ -d "$APP_DIR/.git" ]; then
    echo "Updating repository..."
    cd "$APP_DIR"
    git pull
else
    echo "Cloning repository..."
    git clone "$REPO_URL" "$APP_DIR"
fi

cd "$APP_DIR"

# 4. Configure environment
if [ ! -f .env ]; then
    cp .env.example .env
fi

# Set DOMAIN from instance metadata (if configured)
DOMAIN=$(curl -s -H "Metadata-Flavor: Google" \
    "http://metadata.google.internal/computeMetadata/v1/instance/attributes/domain" 2>/dev/null || echo "")
if [ -n "$DOMAIN" ]; then
    sed -i "s/^DOMAIN=.*/DOMAIN=$DOMAIN/" .env
fi

# 5. Override docker-compose volumes to use persistent disk
export COMPOSE_FILE="docker-compose.yml:docker-compose.gcp.yml"

# 6. Start services
echo "Starting Docker Compose services..."
docker compose pull
docker compose up -d

echo "=== Startup script completed at $(date) ==="
```

### 2. Shutdown Script

The shutdown script handles graceful termination during preemption.

**Metadata key:** `shutdown-script`

**Critical constraint:** GCP gives only 30 seconds for preemption shutdown. The existing worker has `stop_grace_period: 300s` which cannot be honored during preemption.

**Script responsibilities:**
1. Signal Docker Compose to stop gracefully
2. Wait for worker to complete current task (best effort)
3. Ensure filesystem sync

**Example shutdown script:**

```bash
#!/bin/bash

# Log shutdown
echo "=== Shutdown script started at $(date) ===" >> /var/log/shutdown-script.log

# Check if this is a preemption
IS_PREEMPTED=$(curl -s -H "Metadata-Flavor: Google" \
    "http://metadata.google.internal/computeMetadata/v1/instance/preempted" 2>/dev/null || echo "FALSE")

echo "Preemption status: $IS_PREEMPTED" >> /var/log/shutdown-script.log

APP_DIR="/opt/qm-nmr-calc"
cd "$APP_DIR" 2>/dev/null || exit 0

# Graceful shutdown with timeout
# GCP gives ~30s, Docker needs time to stop
export COMPOSE_FILE="docker-compose.yml:docker-compose.gcp.yml"
timeout 25 docker compose down || docker compose kill

# Sync filesystems
sync

echo "=== Shutdown script completed at $(date) ===" >> /var/log/shutdown-script.log
```

### 3. Persistent Disk Configuration

A separate persistent disk stores all Docker volumes, ensuring data survives preemption.

**Disk specifications:**

| Property | Recommended Value | Rationale |
|----------|-------------------|-----------|
| Type | pd-ssd | Fast I/O for NWChem scratch |
| Size | 100 GB minimum | Job data, certificates, logs |
| Zone | Same as VM | Required for attachment |

**Mount strategy:**

```
Persistent Disk (/dev/sdb)
    |
    +-- /mnt/data/
          |
          +-- docker-volumes/
          |     +-- app-data/     -> Docker volume: qm-nmr-calc-data
          |     +-- caddy-data/   -> Docker volume: qm-nmr-calc-caddy-data
          |     +-- caddy-config/ -> Docker volume: qm-nmr-calc-caddy-config
          |
          +-- logs/               -> Application logs
```

### 4. Docker Compose Override (NEW FILE)

Create `docker-compose.gcp.yml` to override volume locations:

```yaml
# docker-compose.gcp.yml
# GCP-specific overrides for persistent disk volumes
# Usage: docker compose -f docker-compose.yml -f docker-compose.gcp.yml up -d

services:
  # Reduce worker stop_grace_period for GCP preemption (30s limit)
  worker:
    stop_grace_period: 25s

volumes:
  app-data:
    driver: local
    driver_opts:
      type: none
      o: bind
      device: /mnt/data/docker-volumes/app-data

  caddy_data:
    driver: local
    driver_opts:
      type: none
      o: bind
      device: /mnt/data/docker-volumes/caddy-data

  caddy_config:
    driver: local
    driver_opts:
      type: none
      o: bind
      device: /mnt/data/docker-volumes/caddy-config
```

---

## Data Persistence Strategy

### What Persists Across Preemption

| Data | Storage | Survives Preemption? |
|------|---------|---------------------|
| Job data (inputs, outputs) | Persistent disk | YES |
| Huey queue state | Persistent disk | YES |
| Let's Encrypt certificates | Persistent disk | YES |
| Docker images | Boot disk cache | NO - repulled |
| Container logs | Ephemeral | NO - lost |
| In-progress calculations | Memory | NO - lost |

### Handling In-Progress Jobs

**Problem:** NWChem calculations can run for 5-30 minutes. A preemption during calculation loses all progress.

**Mitigation strategies:**

1. **Job checkpointing (NOT recommended for NWChem)**
   - NWChem doesn't natively support checkpointing DFT calculations
   - Would require significant application changes

2. **Job retry on restart (RECOMMENDED)**
   - Mark in-progress jobs as "interrupted" on shutdown
   - Re-queue interrupted jobs on startup
   - Requires minor application code change

3. **Accept occasional job loss (SIMPLEST)**
   - User resubmits failed jobs
   - Acceptable for research/non-critical workloads
   - No code changes required

**Recommendation:** Start with option 3 (accept occasional loss), add option 2 (job retry) if preemption frequency is problematic.

### Certificate Persistence

Let's Encrypt certificates are stored in the `caddy_data` volume. With persistent disk storage:

- Certificates survive preemption
- No rate limit issues from repeated certificate requests
- Domain must still point to new IP after preemption (if using ephemeral IP)

**DNS consideration:** Use a static IP or implement dynamic DNS to maintain domain resolution after VM recreation.

---

## VM Configuration

### Recommended Machine Type

For NWChem DFT calculations requiring 32+ cores and 64+ GB RAM:

| Option | vCPUs | RAM | Spot Price (est.) | Notes |
|--------|-------|-----|-------------------|-------|
| n2-highmem-32 | 32 | 256 GB | ~$0.23/hr | Excess RAM, good headroom |
| n2-standard-32 | 32 | 128 GB | ~$0.18/hr | Balanced option |
| n2-custom-32-65536 | 32 | 64 GB | ~$0.15/hr | Exact requirement, custom |
| c2-standard-30 | 30 | 120 GB | ~$0.16/hr | Compute-optimized |

**Recommendation:** `n2-standard-32` provides good balance. Use custom machine type if cost optimization is critical.

### Instance Creation Command

```bash
gcloud compute instances create qm-nmr-calc \
    --zone=us-central1-a \
    --machine-type=n2-standard-32 \
    --provisioning-model=SPOT \
    --instance-termination-action=STOP \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --boot-disk-size=50GB \
    --boot-disk-type=pd-ssd \
    --disk=name=qm-nmr-calc-data,device-name=data-disk,mode=rw,boot=no \
    --metadata-from-file=startup-script=startup.sh,shutdown-script=shutdown.sh \
    --metadata=domain=nmr.example.com \
    --tags=http-server,https-server \
    --scopes=storage-ro
```

### Firewall Rules

```bash
# Allow HTTP/HTTPS (for Caddy)
gcloud compute firewall-rules create allow-http \
    --direction=INGRESS \
    --priority=1000 \
    --network=default \
    --action=ALLOW \
    --rules=tcp:80,tcp:443 \
    --target-tags=http-server,https-server
```

---

## Preemption Handling

### 30-Second Constraint

GCP provides approximately 30 seconds between preemption notice and forced termination.

**Timeline:**

```
T+0s   : Preemption notice (ACPI G2 Soft Off)
         - shutdown-script starts
         - docker compose down issued (SIGINT to containers)

T+0-25s: Graceful shutdown window
         - Worker receives SIGINT
         - Huey attempts to complete current task
         - If task completes, clean exit
         - If task exceeds 25s, interrupted

T+25s  : docker compose down times out, issues SIGKILL
         - Filesystem sync

T+30s  : VM forcibly stopped (ACPI G3 Mechanical Off)
```

### Worker Considerations

The existing worker configuration:

```yaml
worker:
  stop_signal: SIGINT      # Huey graceful shutdown
  stop_grace_period: 300s  # 5 minutes for long calculations
```

**Problem:** The 300s grace period is incompatible with GCP's 30s preemption limit.

**Solution:** Override in `docker-compose.gcp.yml`:

```yaml
worker:
  stop_grace_period: 25s   # GCP preemption compatible
```

**Trade-off:** Jobs running longer than ~20s at preemption time will be interrupted.

### Detecting Preemption

Application code can detect preemption by polling the metadata server:

```python
import requests

def is_preempted():
    try:
        response = requests.get(
            'http://metadata.google.internal/computeMetadata/v1/instance/preempted',
            headers={'Metadata-Flavor': 'Google'},
            timeout=1
        )
        return response.text.strip() == 'TRUE'
    except:
        return False
```

This is useful for logging but not for preventing data loss (insufficient time).

---

## Build Order

### Suggested Implementation Sequence

**Phase 1: Infrastructure Setup (No code changes)**

1. Create persistent disk in GCP
2. Write startup script (metadata)
3. Write shutdown script (metadata)
4. Create `docker-compose.gcp.yml` override file
5. Test: Create spot VM, verify services start

**Phase 2: DNS and HTTPS**

1. Configure static IP or dynamic DNS
2. Set domain in instance metadata
3. Test: Verify Caddy obtains certificate
4. Test: Verify certificate survives preemption

**Phase 3: Robustness (Optional)**

1. Add job retry logic for interrupted calculations
2. Implement preemption logging/monitoring
3. Consider managed instance group for auto-recreation

### What Changes vs. What Stays

| Component | Change Required |
|-----------|-----------------|
| `docker-compose.yml` | NO CHANGE |
| `Caddyfile` | NO CHANGE |
| `Dockerfile.api` | NO CHANGE |
| `Dockerfile.worker` | NO CHANGE |
| `.env` | Minor (domain config) |
| Application code | NO CHANGE (unless job retry needed) |
| NEW: `docker-compose.gcp.yml` | CREATE |
| NEW: `startup.sh` | CREATE |
| NEW: `shutdown.sh` | CREATE |
| NEW: GCP infrastructure | CREATE (disk, VM, firewall) |

---

## Monitoring and Operations

### Logs

```bash
# Startup script logs (on VM)
cat /var/log/startup-script.log

# Shutdown script logs (on VM)
cat /var/log/shutdown-script.log

# Docker Compose logs
docker compose logs -f

# Serial console (GCP Console)
gcloud compute instances get-serial-port-output qm-nmr-calc
```

### Health Checks

The existing Docker health checks work unchanged:

- API: `curl -f http://localhost:8000/health`
- Worker: `pgrep -f huey_consumer`

Add GCP-level monitoring:

```bash
# Check VM status
gcloud compute instances describe qm-nmr-calc --format='get(status)'

# Check if preempted recently
gcloud compute operations list --filter="targetLink:qm-nmr-calc AND operationType:compute.instances.preempted"
```

### Restart After Preemption

If termination action is STOP (recommended), VM enters TERMINATED state:

```bash
# Manual restart
gcloud compute instances start qm-nmr-calc

# Check if resources available
gcloud compute instances start qm-nmr-calc --zone=us-central1-a
```

**Auto-restart option:** Use a Managed Instance Group (MIG) with size=1 for automatic recreation after preemption.

---

## Cost Estimation

### Spot VM Pricing (Approximate)

| Component | On-Demand | Spot (60-70% discount) |
|-----------|-----------|------------------------|
| n2-standard-32 (us-central1) | ~$1.50/hr | ~$0.45/hr |
| 100GB pd-ssd | ~$17/month | $17/month (no discount) |
| Static IP (if used) | $7.20/month | $7.20/month |
| Egress (1TB/month) | ~$85/month | ~$85/month |

**Estimated monthly cost (24/7 operation):**
- VM (spot): ~$330/month
- Storage: ~$17/month
- Network: Variable
- **Total: ~$350-400/month**

Compare to on-demand: ~$1,100/month

### Cost Optimization

1. **Stop when idle:** Implement auto-stop when no jobs pending
2. **Use preemption-tolerant scheduling:** Run during off-peak hours
3. **Consider committed use discounts:** If running consistently, CUDs may beat spot

---

## Alternative Approaches

### Considered but Not Recommended

| Alternative | Why Not |
|-------------|---------|
| Cloud Run | No support for 300s+ jobs, cold start issues |
| GKE with Spot pods | Overkill for single-worker deployment |
| Standard VM | 2-3x more expensive |
| Preemptible VM | Being deprecated, Spot has better features |

### Future Considerations

1. **Managed Instance Group (MIG):** Auto-recreation after preemption
2. **Instance scheduler:** Stop overnight when not in use
3. **Regional persistent disk:** Higher availability if needed
4. **Cloud Filestore:** If multiple workers needed (ReadWriteMany)

---

## Sources

### Official Documentation (HIGH confidence)

- [Spot VMs Overview](https://docs.cloud.google.com/compute/docs/instances/spot)
- [Create and Use Spot VMs](https://docs.cloud.google.com/compute/docs/instances/create-use-spot)
- [Startup Scripts on Linux VMs](https://docs.cloud.google.com/compute/docs/instances/startup-scripts/linux)
- [Run Shutdown Scripts](https://docs.cloud.google.com/compute/docs/shutdownscript)
- [Persistent Disks](https://docs.cloud.google.com/compute/docs/disks/persistent-disks)
- [Stateful Managed Instance Groups](https://docs.cloud.google.com/compute/docs/instance-groups/stateful-migs)
- [Machine Images](https://docs.cloud.google.com/compute/docs/machine-images)

### Supporting Resources (MEDIUM confidence)

- [Run Docker on GCP Compute Instance VMs](https://www.pascallandau.com/blog/gcp-compute-instance-vm-docker/)
- [GCP Spot VMs Explained](https://www.pump.co/blog/spot-instances-gcp)
- [Pro Tip: Use Shutdown Script to Detect Preemption](https://blog.ceshine.net/post/gcp-shutdown-script/)
- [Google Cloud Compute Engine Snapshot Use](https://bluexp.netapp.com/blog/gcp-google-cloud-compute-engine-snapshot-use-cvo-blg)

### Pricing References

- [CloudPrice GCP Compute](https://cloudprice.net/gcp/compute)
- [n2-highmem-32 Pricing](https://www.economize.cloud/resources/gcp/pricing/compute-engine/n2-highmem-32/)
- [GCP Instances Comparison](https://gcpinstances.doit.com/)
