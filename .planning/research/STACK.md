# Stack Research: GCP Spot VM Deployment

**Project:** qm-nmr-calc GCP Spot Instance Deployment
**Researched:** 2026-02-04
**Overall Confidence:** HIGH (verified via official GCP documentation)

## Executive Summary

Deploying the existing Docker Compose stack to GCP Spot VMs requires:
1. **gcloud CLI** for VM provisioning (no Terraform needed for single-VM manual lifecycle)
2. **c2d-highmem-8** machine type (8 vCPU, 64 GB RAM) for NWChem memory requirements
3. **Ubuntu 22.04 LTS** over Container-Optimized OS (simpler Docker Compose support)
4. **Static external IP** with DNS update on boot for Caddy HTTPS
5. **Startup script** to pull and run docker compose on VM creation

Key insight: GCP Spot VMs provide ~80% discount (~$0.10/hour vs ~$0.49/hour for c2d-highmem-8) but can be preempted with only 30 seconds warning. The existing worker's graceful shutdown (SIGINT, 5-minute grace) cannot fully protect against preemption, so job checkpointing should be considered for long calculations.

## Recommended Stack

### Infrastructure Management

| Component | Choice | Version | Rationale |
|-----------|--------|---------|-----------|
| Provisioning | gcloud CLI | >=555.0.0 | Simple single-VM deployment; Terraform overkill for manual start/stop lifecycle |
| Authentication | gcloud auth | (built-in) | Application Default Credentials for local scripts |
| Project Setup | gcloud projects | (built-in) | One-time project/billing setup |

**Why gcloud CLI over Terraform:**
- Single VM with manual lifecycle (start/stop by user)
- No auto-scaling, no complex state management needed
- Faster iteration during development
- Can always migrate to Terraform later if needed

### VM Instance Configuration

| Component | Choice | Rationale |
|-----------|--------|-----------|
| Machine Type | c2d-highmem-8 | 8 vCPU, 64 GB RAM; NWChem needs ~2GB/MPI process, allows NWCHEM_NPROC=8 with headroom |
| Provisioning Model | Spot | ~80% discount; acceptable for batch workloads with manual lifecycle |
| Termination Action | STOP | Preserves disk on preemption; can restart manually |
| Boot Disk | Ubuntu 22.04 LTS, 50GB SSD | Native Docker Compose support; sufficient for images + job data |
| Zone | us-central1-a | Good spot availability; adjust based on user location |

**Why c2d-highmem-8 over c2-standard-8:**
- c2d-highmem-8: 8 vCPU, **64 GB** RAM (~$0.097/hour spot)
- c2-standard-8: 8 vCPU, **32 GB** RAM (~$0.088/hour spot)
- NWChem recommends 2GB/MPI process minimum; 64GB allows comfortable 8-process parallelism with memory for conformer search and overhead

**Why Ubuntu over Container-Optimized OS:**
- Ubuntu has native `docker compose` support (apt install docker-compose-plugin)
- COS requires workarounds (docker-in-docker, credential helpers)
- Ubuntu allows package installation for debugging if needed
- COS benefits (security hardening, auto-updates) less relevant for manual lifecycle VM

### Networking

| Component | Choice | Rationale |
|-----------|--------|-----------|
| External IP | Static (reserved) | Stable IP for DNS A record; survives stop/start |
| Firewall Rules | HTTP (80), HTTPS (443), SSH (22) | Caddy auto-HTTPS, secure access |
| DNS | External provider or Cloud DNS | A record pointing to static IP |

**Static IP rationale:**
- Spot VMs with STOP termination action retain their IP on restart
- Static IP ($7.30/month when attached) cheaper than DNS propagation delays
- Caddy obtains certificates on first boot; static IP ensures domain resolution

### Cost Estimates (us-central1)

| Machine Type | vCPU | RAM | On-Demand | Spot | Monthly (Spot, 8hr/day) |
|--------------|------|-----|-----------|------|-------------------------|
| c2-standard-8 | 8 | 32 GB | $0.40/hr | $0.088/hr | ~$21 |
| c2d-highmem-8 | 8 | 64 GB | $0.49/hr | $0.097/hr | ~$23 |
| c2-standard-16 | 16 | 64 GB | $0.80/hr | $0.177/hr | ~$42 |
| c2d-highmem-16 | 16 | 128 GB | $0.98/hr | $0.195/hr | ~$47 |

**Recommendation:** c2d-highmem-8 at ~$0.10/hour spot. At 8 hours/day typical use: ~$23/month. With static IP: ~$30/month total.

**Cost comparison with always-on VPS:**
- DigitalOcean 8 vCPU, 16GB: $96/month (always-on)
- GCP c2d-highmem-8 spot, 8hr/day: ~$30/month (when needed)

## Deployment Architecture

### Startup Script Flow

```
VM Start
    |
    v
[cloud-init / startup-script]
    |
    +-- Install Docker (if not present)
    +-- Install Docker Compose plugin
    +-- Clone/pull qm-nmr-calc repo
    +-- Copy .env from GCS or metadata
    +-- docker compose pull
    +-- docker compose up -d
    |
    v
[Caddy obtains HTTPS certificate]
    |
    v
[Service ready]
```

### Shutdown Script Flow (Preemption Handling)

```
Preemption Notice (30 seconds)
    |
    v
[shutdown-script]
    |
    +-- docker compose down (SIGINT to worker)
    +-- Worker has 30s to finish current task
    |   (Note: 5-min grace period truncated by preemption)
    +-- Optional: Upload job state to GCS
    |
    v
[VM Stopped]
```

**Important limitation:** The worker's 5-minute stop_grace_period cannot be honored during preemption. Only ~25 seconds available after shutdown script starts. Long-running NWChem calculations may be interrupted.

### Manual Lifecycle Commands

```bash
# Start (create or restart stopped VM)
gcloud compute instances start qm-nmr-calc --zone=us-central1-a

# Stop (graceful, preserves disk)
gcloud compute instances stop qm-nmr-calc --zone=us-central1-a

# Check status
gcloud compute instances describe qm-nmr-calc --zone=us-central1-a \
    --format="value(status)"

# SSH for debugging
gcloud compute ssh qm-nmr-calc --zone=us-central1-a

# View logs
gcloud compute ssh qm-nmr-calc --zone=us-central1-a \
    --command="docker compose logs -f"
```

## gcloud CLI Commands Reference

### One-Time Setup

```bash
# Install gcloud CLI (macOS)
brew install google-cloud-sdk

# Or download installer
curl https://sdk.cloud.google.com | bash

# Initialize and authenticate
gcloud init
gcloud auth login

# Set default project
gcloud config set project YOUR_PROJECT_ID

# Enable required APIs
gcloud services enable compute.googleapis.com
```

### Reserve Static IP

```bash
gcloud compute addresses create qm-nmr-calc-ip \
    --region=us-central1 \
    --description="Static IP for qm-nmr-calc"

# Get the IP address
gcloud compute addresses describe qm-nmr-calc-ip \
    --region=us-central1 \
    --format="value(address)"
```

### Create Firewall Rules

```bash
# Allow HTTP/HTTPS for Caddy
gcloud compute firewall-rules create allow-http-https \
    --direction=INGRESS \
    --action=ALLOW \
    --rules=tcp:80,tcp:443 \
    --source-ranges=0.0.0.0/0 \
    --target-tags=http-server

# Allow SSH (consider restricting source-ranges)
gcloud compute firewall-rules create allow-ssh \
    --direction=INGRESS \
    --action=ALLOW \
    --rules=tcp:22 \
    --source-ranges=0.0.0.0/0 \
    --target-tags=ssh-server
```

### Create Spot VM

```bash
gcloud compute instances create qm-nmr-calc \
    --zone=us-central1-a \
    --machine-type=c2d-highmem-8 \
    --provisioning-model=SPOT \
    --instance-termination-action=STOP \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --boot-disk-size=50GB \
    --boot-disk-type=pd-ssd \
    --address=qm-nmr-calc-ip \
    --tags=http-server,ssh-server \
    --metadata-from-file=startup-script=startup.sh \
    --scopes=storage-ro
```

### Startup Script Template

```bash
#!/bin/bash
set -e

# Install Docker
curl -fsSL https://get.docker.com | sh

# Install Docker Compose plugin
apt-get update
apt-get install -y docker-compose-plugin

# Clone or update repository
if [ -d "/opt/qm-nmr-calc" ]; then
    cd /opt/qm-nmr-calc && git pull
else
    git clone https://github.com/steinbeck/qm-nmr-calc.git /opt/qm-nmr-calc
fi

cd /opt/qm-nmr-calc

# Configure environment (from GCS or hardcoded)
cat > .env << 'EOF'
DOMAIN=nmr.example.com
ACME_EMAIL=admin@example.com
NWCHEM_NPROC=8
OMP_NUM_THREADS=8
WORKER_MEMORY_LIMIT=48g
EOF

# Pull latest images and start
docker compose pull
docker compose up -d

# Log startup completion
echo "qm-nmr-calc started at $(date)" >> /var/log/qm-nmr-calc-startup.log
```

## Alternatives Considered

### Infrastructure as Code

| Option | Why Not Chosen |
|--------|----------------|
| Terraform | Overkill for single VM with manual lifecycle; adds state management complexity |
| Pulumi | Same as Terraform; better for multi-resource, versioned infrastructure |
| Cloud Deployment Manager | GCP-specific, less portable; gcloud CLI simpler for this use case |
| Ansible | Good for configuration, but gcloud handles VM creation natively |

**When to reconsider Terraform:**
- If adding auto-scaling or managed instance groups
- If deploying to multiple environments (dev/staging/prod)
- If other team members need to manage infrastructure

### Operating System

| Option | Why Not Chosen |
|--------|----------------|
| Container-Optimized OS | Docker Compose requires workarounds; no apt for debugging |
| Debian | Ubuntu has better Docker Compose documentation and support |
| Fedora CoreOS | Designed for Kubernetes; overkill for Docker Compose |

### VM Provisioning Model

| Option | Why Not Chosen |
|--------|----------------|
| On-Demand | 5x more expensive; not justified for batch workloads |
| Committed Use | Requires 1-3 year commitment; user wants on-demand start/stop |
| Preemptible (legacy) | Spot VMs are the newer, recommended replacement |

### Machine Type Family

| Option | Why Not Chosen |
|--------|----------------|
| N2/N2D (General) | C2/C2D optimized for compute-heavy NWChem workloads |
| E2 (Cost-optimized) | Shared-core options; NWChem needs dedicated CPU |
| M1/M2 (Memory) | Overkill; 64GB from c2d-highmem sufficient |

## Pitfalls and Mitigations

### Preemption During Long Calculations

**Risk:** NWChem calculations can take 30+ minutes. Preemption gives only 30 seconds.

**Mitigations:**
1. Use STOP termination action (preserves disk state)
2. Implement job checkpointing in future milestone
3. Consider on-demand for critical long-running jobs
4. Monitor preemption rates in chosen zone

### Certificate Acquisition on First Boot

**Risk:** Caddy needs domain to resolve to VM IP before obtaining Let's Encrypt cert.

**Mitigations:**
1. Use static IP reserved before VM creation
2. Configure DNS before first boot
3. Set low TTL (300s) on A record for quick updates
4. First boot may take 1-2 minutes for certificate acquisition

### Docker Image Pull Times

**Risk:** First boot pulls ~2GB of worker images, adding startup latency.

**Mitigation:**
1. Pre-pull images on first setup, then use `docker compose pull` for updates
2. Consider custom VM image with pre-installed images (advanced)
3. Use boot disk snapshot after first successful setup

### Cost Accumulation When Stopped

**Risk:** Static IP charges accumulate when VM stopped ($7.30/month).

**Mitigation:**
1. Accept as cost of stable DNS (cheaper than alternatives)
2. Or use dynamic DNS with ephemeral IP (more complex)

## Sources

### Official Documentation (HIGH confidence)
- [GCP Spot VMs](https://docs.cloud.google.com/compute/docs/instances/spot)
- [Create and Use Spot VMs](https://docs.cloud.google.com/compute/docs/instances/create-use-spot)
- [Compute-Optimized Machines](https://docs.cloud.google.com/compute/docs/compute-optimized-machines)
- [gcloud compute instances create](https://cloud.google.com/sdk/gcloud/reference/compute/instances/create)
- [Configure Static External IP](https://docs.cloud.google.com/compute/docs/ip-addresses/configure-static-external-ip-address)
- [VPC Firewall Rules](https://docs.cloud.google.com/firewall/docs/using-firewalls)
- [Graceful Shutdown Overview](https://docs.cloud.google.com/compute/docs/instances/graceful-shutdown-overview)
- [gcloud CLI Release Notes](https://docs.cloud.google.com/sdk/docs/release-notes) - Version 555.0.0

### Pricing (MEDIUM confidence - prices vary by region and time)
- [Vantage c2-standard-8](https://instances.vantage.sh/gcp/c2-standard-8) - $0.088/hr spot
- [CloudPrice c2d-highmem-8](https://cloudprice.net/gcp/compute/instances/c2d-highmem-8) - $0.097/hr spot
- [GCP VM Instance Pricing](https://cloud.google.com/compute/vm-instance-pricing)

### Community Patterns (MEDIUM confidence)
- [Docker on GCP Compute VMs](https://www.pascallandau.com/blog/gcp-compute-instance-vm-docker/)
- [Docker Compose on Container-Optimized OS](https://gist.github.com/robvanoostenrijk/d14e67c0b89e65a429d05bf4086d1947)
- [GCP Shutdown Script for Preemption](https://blog.ceshine.net/post/gcp-shutdown-script/)
- [Stopping Docker Container on COS](https://medium.com/google-cloud/stopping-a-docker-contain-on-cos-9a1f615dc85)

## Roadmap Implications

### Recommended Phase Structure

1. **Phase 1: GCP Setup** - Project, billing, static IP, firewall rules, DNS
2. **Phase 2: VM Deployment Script** - gcloud commands, startup script, .env management
3. **Phase 3: Lifecycle Scripts** - Start/stop wrapper scripts, status checking
4. **Phase 4: Preemption Handling** - Shutdown script, job state preservation (optional)

### Research Flags for Implementation

- **DNS Configuration:** User needs to configure their domain provider; document common providers
- **Secrets Management:** .env file handling - GCS bucket, metadata, or hardcoded in startup script
- **Job Recovery:** If preemption during calculation is a concern, may need deeper research on NWChem restart files
