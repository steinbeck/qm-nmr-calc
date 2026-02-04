# Research Summary: v2.6 GCP Spot VM Deployment

**Project:** qm-nmr-calc v2.6 - Google Cloud Spot Deployment
**Domain:** Docker application deployment to GCP Spot VMs
**Researched:** 2026-02-04
**Confidence:** HIGH

## Executive Summary

GCP Spot VM deployment for qm-nmr-calc is a well-understood pattern with one critical constraint: **GCP provides only 30 seconds warning before preemption**, which is incompatible with NWChem calculations that run for hours. The recommended approach accepts this limitation for v2.6 (manual lifecycle) and focuses on data persistence and clean shutdown rather than calculation preservation. Research indicates this is the standard pattern for long-running batch workloads on spot infrastructure.

The deployment architecture is straightforward: the existing Docker Compose stack runs unchanged, wrapped by GCP-specific startup/shutdown scripts and persistent disk storage. The gcloud CLI is sufficient for single-VM manual lifecycle management - Terraform adds complexity without benefit for this use case. Ubuntu 22.04 LTS is preferred over Container-Optimized OS for simpler Docker Compose support.

Key risk is cost surprise from static IP charges when idle ($7.30/month) and potential calculation loss during preemption. Both are acceptable trade-offs: the static IP enables stable DNS for Caddy HTTPS, and preemption-interrupted jobs can simply be resubmitted. At ~$0.10/hour spot pricing vs ~$0.49/hour on-demand (c2d-highmem-8), the 80% discount justifies accepting occasional interruption.

## Key Findings

### Recommended Stack

The gcloud CLI provides everything needed for single-VM deployment without the state management overhead of Terraform. Key configuration choices are driven by NWChem's resource requirements and Caddy's HTTPS needs.

**Core technologies:**
- **gcloud CLI (>=555.0.0):** VM provisioning, lifecycle management - simpler than Terraform for manual start/stop
- **c2d-highmem-8:** 8 vCPU, 64 GB RAM - NWChem needs ~2GB/MPI process; allows 8-process parallelism with headroom
- **Ubuntu 22.04 LTS:** Native docker compose plugin support; COS requires workarounds
- **Persistent disk (pd-ssd, 50-100GB):** Job data survives preemption; separate from boot disk
- **Static external IP:** Stable DNS for Caddy HTTPS; survives stop/start

**Cost estimate:** ~$30/month at 8 hours/day usage (c2d-highmem-8 spot + static IP)

### Expected Features

**Must have (table stakes):**
- One-command VM creation with all configuration
- Spot VM with STOP termination action (preserves disk)
- Startup script: Docker install, repo clone, docker compose up
- Shutdown script: graceful container stop (25s limit)
- Firewall rules for HTTP/HTTPS/SSH

**Should have (differentiators):**
- Lifecycle scripts: stop, start, status, delete, ssh, logs
- Interactive deployment with defaults (zone, machine type)
- Configuration persistence (.gcp-instance file)
- Cost estimate display during deployment

**Defer (v2+):**
- Terraform/Pulumi IaC
- Auto-scaling / managed instance groups
- Cloud DNS integration
- Job checkpointing for preemption recovery
- Monitoring dashboards

### Architecture Approach

The Docker Compose stack remains **completely unchanged**. GCP deployment wraps it with infrastructure: startup script handles Docker installation and service startup on boot; shutdown script triggers graceful container shutdown during preemption; persistent disk stores all Docker volumes. A new `docker-compose.gcp.yml` override file reduces worker stop_grace_period from 300s to 25s for GCP compatibility.

**Major components:**
1. **Deployment script (deploy-gcp.sh):** Creates VM, firewall rules, persistent disk via gcloud CLI
2. **Startup script (startup.sh):** Runs on every boot - mounts disk, installs Docker, starts services
3. **Shutdown script (shutdown.sh):** Triggered on preemption - graceful docker compose down within 25s
4. **docker-compose.gcp.yml:** Volume overrides for persistent disk, reduced stop_grace_period
5. **Lifecycle scripts:** stop-gcp.sh, start-gcp.sh, status-gcp.sh, delete-gcp.sh

### Critical Pitfalls

1. **30-second preemption window vs hours-long calculations** - Accept job loss for v2.6; document clearly; users resubmit. Future: implement job retry logic.

2. **Using local SSD instead of persistent disk** - All data lost on preemption. Always use persistent disk for job data; mount at /mnt/data.

3. **Containers not starting after restart** - Use `restart: unless-stopped` policy; startup script runs `docker compose up -d` idempotently on every boot.

4. **Firewall rules missing in custom VPC** - Use default VPC (has pre-configured rules) or create HTTP/HTTPS/SSH rules explicitly before VM creation.

5. **DNS pointing to old IP after restart** - Use reserved static IP (not ephemeral); IP survives stop/start. $7.30/month cost is acceptable.

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: GCP Infrastructure Setup
**Rationale:** Foundation for all other phases; firewall rules and disk must exist before VM
**Delivers:** GCP project configured, static IP reserved, firewall rules created, persistent disk created
**Addresses:** Network access (table stake), security basics (table stake)
**Avoids:** Pitfall #4 (firewall rules missing)

### Phase 2: VM Deployment Script
**Rationale:** Core deliverable of v2.6; depends on infrastructure from Phase 1
**Delivers:** deploy-gcp.sh that creates Spot VM with all configuration
**Uses:** gcloud CLI, startup script, shutdown script
**Implements:** One-command deployment (table stake)
**Avoids:** Pitfall #6 (--preemptible flag), Pitfall #7 (DELETE termination action)

### Phase 3: Lifecycle Management Scripts
**Rationale:** Users need to start/stop/check status after initial deployment
**Delivers:** stop-gcp.sh, start-gcp.sh, status-gcp.sh, delete-gcp.sh, ssh-gcp.sh, logs-gcp.sh
**Addresses:** Stop/start commands (differentiator), operational convenience
**Avoids:** Manual gcloud command memorization

### Phase 4: Documentation and Testing
**Rationale:** Users need clear guidance on prerequisites, costs, limitations
**Delivers:** Deployment guide, cost estimates, troubleshooting guide
**Addresses:** Preemption behavior explanation, machine type recommendations
**Avoids:** Pitfall #1 (surprise at calculation loss during preemption)

### Phase Ordering Rationale

- **Phase 1 before Phase 2:** VM creation fails without firewall rules and disk
- **Phase 2 before Phase 3:** Lifecycle scripts assume VM exists
- **Phase 3 before Phase 4:** Need working deployment to document
- Each phase produces testable artifact before next begins

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 2 (deployment script):** May need research on cloud-init vs metadata startup scripts for reliability

Phases with standard patterns (skip research-phase):
- **Phase 1 (infrastructure):** Well-documented gcloud commands, no complexity
- **Phase 3 (lifecycle scripts):** Simple gcloud wrappers, standard pattern
- **Phase 4 (documentation):** No research needed, document what was built

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Official GCP documentation, verified pricing sources |
| Features | HIGH | Standard deployment patterns, clear user needs |
| Architecture | HIGH | Existing Docker stack unchanged; GCP patterns well-documented |
| Pitfalls | HIGH | Official docs + community best practices confirm issues |

**Overall confidence:** HIGH

### Gaps to Address

- **DNS provider integration:** Documentation should cover common providers (Cloudflare, Namecheap, Google Domains) but not automate - user handles DNS manually
- **Secrets management:** .env file handling options (GCS bucket, metadata, hardcoded) need decision during implementation
- **NWChem restart files:** If preemption-during-calculation becomes a pain point, may need deeper research on permanent_dir configuration for future milestone

## Sources

### Primary (HIGH confidence)
- [GCP Spot VMs Documentation](https://docs.cloud.google.com/compute/docs/instances/spot) - Preemption behavior, pricing
- [Create and Use Spot VMs](https://docs.cloud.google.com/compute/docs/instances/create-use-spot) - gcloud CLI commands
- [Compute-Optimized Machines](https://docs.cloud.google.com/compute/docs/compute-optimized-machines) - Machine type selection
- [GCP Startup Scripts](https://cloud.google.com/compute/docs/instances/startup-scripts/linux) - Boot automation
- [GCP Shutdown Scripts](https://cloud.google.com/compute/docs/shutdownscript) - Graceful preemption handling
- [Persistent Disk Documentation](https://cloud.google.com/compute/docs/disks/persistent-disks) - Data persistence

### Secondary (MEDIUM confidence)
- [CloudPrice c2d-highmem-8](https://cloudprice.net/gcp/compute/instances/c2d-highmem-8) - Spot pricing estimates
- [Docker on GCP Compute VMs](https://www.pascallandau.com/blog/gcp-compute-instance-vm-docker/) - Deployment patterns
- [GCP Shutdown Script Blog](https://blog.ceshine.net/post/gcp-shutdown-script/) - Preemption detection

### Tertiary (LOW confidence)
- NWChem restart file behavior during preemption - needs validation if job checkpointing pursued in future

---
*Research completed: 2026-02-04*
*Ready for roadmap: yes*
