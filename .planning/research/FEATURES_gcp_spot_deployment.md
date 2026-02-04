# Feature Landscape: GCP Spot VM Deployment

**Domain:** Docker app deployment to GCP Spot VMs
**Researched:** 2026-02-04
**Confidence:** HIGH (GCP official documentation, established patterns)
**Target Users:** Researchers who need cost-effective cloud compute for expensive NMR calculations

## Executive Summary

GCP Spot VM deployment for qm-nmr-calc addresses a specific use case: **on-demand, cost-optimized cloud compute** for researchers who need more CPU/memory than their local machines but don't want to pay for always-on servers.

The existing Docker Compose stack (worker, API, Caddy) is well-suited for Spot VM deployment because:
1. **Graceful shutdown handling already exists** (SIGINT, 5-minute grace period)
2. **Health checks are in place** for service recovery
3. **Stateless calculations** - job data persists in volumes, no distributed state to manage

The core feature set focuses on **one-command deployment** that handles VM provisioning, Docker setup, and service startup. Since this is manual lifecycle (user starts/stops as needed), the deployment tooling is simpler than always-on infrastructure.

## Table Stakes Features

Features users **expect** from cloud deployment. Missing these = deployment feels broken or incomplete.

### VM Provisioning

| Feature | Why Expected | Complexity | Depends On | Notes |
|---------|--------------|------------|------------|-------|
| **One-command VM creation** | Core promise of the milestone | MEDIUM | gcloud CLI | `gcloud compute instances create` with all options |
| **Spot VM configuration** | 60-91% cost savings vs on-demand | LOW | gcloud CLI | `--provisioning-model=SPOT` flag |
| **Appropriate machine type** | NWChem needs high CPU/memory | LOW | gcloud CLI | c2-standard-8 or n2-highmem-4 recommended |
| **Boot disk sizing** | Docker images are ~2-3 GB | LOW | gcloud CLI | 40 GB minimum for worker image |
| **Network access** | Web UI must be reachable | LOW | Firewall rules | HTTP/HTTPS from 0.0.0.0/0 |

**Rationale:** Users expect a single script/command that provisions a working VM. Manual console clicking or multi-step processes defeat the purpose of automation.

### Docker Setup (Startup Script)

| Feature | Why Expected | Complexity | Depends On | Notes |
|---------|--------------|------------|------------|-------|
| **Docker installation** | Required for containers | LOW | Startup script | `curl -fsSL https://get.docker.com | sh` |
| **Docker Compose installation** | Required for multi-container stack | LOW | Docker | Included with Docker Engine |
| **Repository clone** | Get docker-compose.yml and Caddyfile | LOW | Git | Clone from GitHub |
| **Service startup** | Containers must run after boot | LOW | Docker Compose | `docker compose up -d` |
| **Environment configuration** | Domain, resource settings | LOW | .env file | Write .env before compose up |

**Rationale:** Startup script must handle complete setup from bare VM to running service. User should not SSH in to run commands manually.

### Preemption Handling

| Feature | Why Expected | Complexity | Depends On | Notes |
|---------|--------------|------------|------------|-------|
| **Graceful container shutdown** | Don't lose in-progress calculations | MEDIUM | cloud-init | SIGINT to containers before VM stops |
| **30-second warning usage** | Maximum time to finish graceful stop | LOW | GCP preemption signal | Already have 5-min grace period, but only 30s available |
| **Job state preservation** | Resume after restart | LOW | Volume persistence | Job data survives preemption |

**Rationale:** GCP provides 30-second warning before preemption. The existing graceful shutdown (SIGINT to Huey worker) must be triggered during this window.

### Access and Security

| Feature | Why Expected | Complexity | Depends On | Notes |
|---------|--------------|------------|------------|-------|
| **Firewall rules for HTTP/HTTPS** | Web UI access | LOW | gcloud CLI | --tags=http-server,https-server |
| **SSH access for debugging** | Troubleshoot issues | LOW | Default | GCP default allows SSH |
| **Service account (minimal)** | Security best practice | LOW | gcloud CLI | Default compute service account is fine |

**Rationale:** Security basics - allow only necessary traffic, use default service accounts rather than overprivileged accounts.

## Differentiators

Features that provide **better UX** over minimal deployment. Not strictly required, but significantly improve experience.

### Deployment Experience

| Feature | Value Proposition | Complexity | Depends On | Notes |
|---------|-------------------|------------|------------|-------|
| **Single deploy script** | `./deploy-gcp.sh` does everything | MEDIUM | gcloud CLI | Script orchestrates all gcloud commands |
| **Interactive prompts for options** | User chooses zone, machine type | LOW | Bash | read -p for user input with defaults |
| **Deployment status output** | Know when ready to use | LOW | Bash | Print URL when deployment completes |
| **Estimated cost display** | Know what you're spending | LOW | Hardcoded estimates | Show hourly cost for chosen machine type |

**Rationale:** Good CLI UX - script guides user through options, shows progress, and provides URL when done.

### Lifecycle Management

| Feature | Value Proposition | Complexity | Depends On | Notes |
|---------|-------------------|------------|------------|-------|
| **Stop command** | Pause VM without deleting | LOW | gcloud CLI | `gcloud compute instances stop` |
| **Start command** | Resume stopped VM | LOW | gcloud CLI | `gcloud compute instances start` |
| **Delete command** | Clean teardown | LOW | gcloud CLI | `gcloud compute instances delete` |
| **Status command** | Check if running | LOW | gcloud CLI | `gcloud compute instances describe` |

**Rationale:** Complete lifecycle: deploy, stop, start, delete. User controls costs by stopping when not in use.

### Operational Convenience

| Feature | Value Proposition | Complexity | Depends On | Notes |
|---------|-------------------|------------|------------|-------|
| **IP address retrieval** | Know how to access | LOW | gcloud CLI | Show external IP after create |
| **SSH shortcut** | Easy debugging access | LOW | gcloud CLI | `gcloud compute ssh [instance]` |
| **Log streaming** | View service logs remotely | LOW | SSH + docker | `gcloud compute ssh ... -- docker compose logs -f` |
| **Configuration persistence** | Remember settings for stop/start | LOW | Local file | Store instance name, zone in .gcp-instance |

**Rationale:** Reduce friction for common operations. User shouldn't need to remember instance names or look up IPs.

### Cost Optimization

| Feature | Value Proposition | Complexity | Depends On | Notes |
|---------|-------------------|------------|------------|-------|
| **Region selection guidance** | Spot prices vary by region | LOW | Documentation | Note which regions have better availability |
| **Instance termination action** | STOP vs DELETE on preemption | LOW | gcloud CLI | `--instance-termination-action=STOP` preserves disk |
| **Automatic restart = false** | Don't auto-restart after preemption | LOW | gcloud CLI | User manually restarts when ready |

**Rationale:** Spot VMs require understanding of preemption behavior. Clear options let user choose appropriate settings.

## Anti-Features

Features to **deliberately NOT build** for this scope. Common overengineering mistakes.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Terraform/Pulumi IaC** | Overkill for single VM, manual lifecycle | Simple gcloud CLI script |
| **Auto-scaling / MIGs** | Manual lifecycle means user controls when running | Single VM only |
| **Load balancer** | Single VM doesn't need it | Caddy handles HTTPS directly |
| **Cloud DNS integration** | Users can point their own domain | Document manual DNS setup |
| **Google-managed SSL certs** | Requires load balancer (cost) | Caddy auto-HTTPS with Let's Encrypt |
| **Persistent disk snapshots** | Job data is regenerable, not precious | User can backup volumes manually |
| **Cloud Monitoring integration** | Overkill for on-demand usage | `docker compose logs` is sufficient |
| **Preemption prediction** | Unreliable and overcomplicated | Accept preemption, handle gracefully |
| **Multi-zone redundancy** | Single VM deployment | User picks one zone |
| **Automatic restart after preemption** | User controls when compute runs | Manual start after preemption |
| **Cloud Run / Cloud Functions** | Long-running calculations don't fit serverless | Stick with VM |
| **GKE / Kubernetes** | Massive overkill for single-service deployment | Docker Compose on single VM |
| **Reserved IP addresses** | Costs money when VM is stopped | Ephemeral IP is fine |
| **Custom images** | Startup script is simpler, more maintainable | Use standard Ubuntu + startup script |

**Rationale:** This is a **manual lifecycle, cost-optimized deployment**. Every additional GCP service adds complexity and cost. Keep it simple: one VM, one script, user controls when it runs.

## Feature Dependencies

```
Prerequisites:
├── gcloud CLI installed locally
├── GCP project with billing enabled
├── Compute Engine API enabled
└── User authenticated (`gcloud auth login`)

Deployment Script (deploy-gcp.sh):
├── VM Creation
│   ├── Machine type selection
│   ├── Zone selection
│   ├── Spot VM configuration
│   └── Network tags (http-server, https-server)
├── Firewall Rules
│   └── Ensure HTTP/HTTPS rules exist (create if missing)
├── Startup Script
│   ├── Docker installation
│   ├── Repo clone
│   ├── .env configuration
│   └── docker compose up -d
└── Output
    ├── External IP
    └── Access URL

cloud-init Configuration:
├── Graceful shutdown handling
│   └── systemd service for docker compose down
└── Container signal forwarding

Lifecycle Scripts:
├── stop-gcp.sh
├── start-gcp.sh
├── status-gcp.sh
└── delete-gcp.sh
```

## Dependencies on Existing Features

The GCP deployment builds directly on existing Docker deployment work:

| Existing Feature | How GCP Deployment Uses It |
|------------------|---------------------------|
| docker-compose.yml | Deployed via startup script, unchanged |
| Caddyfile | Handles HTTPS, unchanged |
| GHCR images | Pulled on VM startup |
| Graceful shutdown (SIGINT) | Triggered during preemption |
| Health checks | Used by Docker restart policies |
| .env configuration | Generated by deploy script |
| 5-minute stop_grace_period | Still relevant, though preemption only gives 30s |

**No changes to existing Docker stack required.** GCP deployment is pure addition of deployment tooling.

## MVP Recommendation

### Must Have (Table Stakes)

1. **deploy-gcp.sh** - Single script that:
   - Prompts for zone, machine type (with defaults)
   - Creates Spot VM with appropriate configuration
   - Runs startup script for Docker setup
   - Outputs access URL when complete

2. **Startup script** - Embedded in deploy script or separate file:
   - Install Docker
   - Clone repo
   - Configure .env with domain/IP
   - Start containers

3. **Firewall rule creation** - Script creates rules if missing:
   - Allow HTTP (80/tcp) from anywhere
   - Allow HTTPS (443/tcp) from anywhere

4. **Lifecycle commands**:
   - `./gcp-stop.sh` - Stop VM (preserves disk)
   - `./gcp-start.sh` - Start stopped VM
   - `./gcp-delete.sh` - Delete VM completely
   - `./gcp-status.sh` - Check VM status

5. **Documentation**:
   - Prerequisites (gcloud setup)
   - Machine type recommendations
   - Cost estimates
   - Preemption behavior explanation

### Should Have (Differentiators)

1. **cloud-init for graceful shutdown** - Proper SIGINT to containers on preemption
2. **Cost display** - Show estimated hourly cost during deployment
3. **SSH shortcut** - Include command to SSH into instance
4. **Log streaming** - Command to view container logs remotely

### Defer to Post-MVP

- Terraform configuration
- Multi-zone support
- Automatic domain configuration
- Cloud DNS integration
- Monitoring dashboards
- Backup automation

## Machine Type Recommendations

Based on NWChem computational requirements:

| Use Case | Machine Type | vCPUs | Memory | Spot Price/hr* | Notes |
|----------|--------------|-------|--------|----------------|-------|
| Light testing | n2-standard-4 | 4 | 16 GB | ~$0.03-0.05 | Quick tests, small molecules |
| Standard use | c2-standard-8 | 8 | 32 GB | ~$0.08-0.15 | Most NMR calculations |
| Heavy use | c2-standard-16 | 16 | 64 GB | ~$0.15-0.30 | Large molecules, conformer ensembles |
| Memory-intensive | n2-highmem-8 | 8 | 64 GB | ~$0.08-0.15 | When memory is bottleneck |

*Spot prices vary by region and time. These are rough estimates for US regions.

**Recommended default:** c2-standard-8 - good balance of CPU (8 cores for MPI) and memory (32 GB).

## Preemption Handling Detail

GCP Spot VM preemption sequence:

```
1. GCP sends ACPI G2 signal (30-second warning)
2. Shutdown scripts run (if configured)
3. After 30 seconds, GCP forces termination (ACPI G3)
```

For qm-nmr-calc:

```
Current graceful shutdown:
- Worker: SIGINT → Huey finishes current task → clean exit
- API: SIGTERM → FastAPI shutdown → clean exit
- Grace period: 5 minutes configured

Problem:
- GCP only provides 30 seconds
- Long NWChem calculations (30+ minutes) cannot complete

Solution:
- Accept that preemption may interrupt calculations
- Job state is preserved (input files exist)
- User can re-submit job after restart
- Document this limitation clearly
```

**Key insight:** Don't try to prevent interruption; handle it gracefully. The goal is clean container shutdown (no corrupt state), not calculation completion.

## Sources

### GCP Spot VMs
- [GCP Spot VM Documentation](https://docs.cloud.google.com/compute/docs/instances/spot) - Official reference, preemption behavior
- [Create and Use Spot VMs](https://cloud.google.com/compute/docs/instances/create-use-spot) - gcloud CLI commands
- [GCP Spot VM Best Practices](https://cloud.google.com/blog/products/compute/google-cloud-spot-vm-use-cases-and-best-practices) - Use cases, cost savings

### Graceful Shutdown
- [Graceful Shutdown Overview](https://docs.cloud.google.com/compute/docs/instances/graceful-shutdown-overview) - 30-second warning, shutdown scripts
- [Shutdown Scripts](https://cloud.google.com/compute/docs/shutdownscript) - How to configure
- [Docker Container Shutdown on COS](https://medium.com/google-cloud/stopping-a-docker-contain-on-cos-9a1f615dc85) - cloud-init approach for Docker

### Startup Scripts and cloud-init
- [Startup Scripts on Linux VMs](https://cloud.google.com/compute/docs/instances/startup-scripts/linux) - Official startup script guide
- [Container-Optimized OS cloud-init](https://cloud.google.com/container-optimized-os/docs/how-to/run-container-instance) - cloud-init for Docker
- [Creating and Configuring COS Instances](https://cloud.google.com/container-optimized-os/docs/how-to/create-configure-instance) - COS configuration

### Firewall and Networking
- [VPC Firewall Rules](https://cloud.google.com/firewall/docs/firewalls) - Network access configuration
- [Firewall Rules for Containers](https://cloud.google.com/compute/docs/containers/configuring-options-to-run-containers) - Port publishing

### Machine Types and Pricing
- [Machine Families Comparison](https://cloud.google.com/compute/docs/machine-resource) - Machine type selection
- [Spot VM Pricing](https://cloud.google.com/spot-vms/pricing) - Current Spot prices by region
- [VM Instance Pricing](https://cloud.google.com/compute/vm-instance-pricing) - On-demand vs Spot comparison

### Deployment Automation
- [Deploy Docker on GCP VMs](https://www.pascallandau.com/blog/gcp-compute-instance-vm-docker/) - Practical tutorial
- [Terraform GCP Spot VMs](https://medium.com/@ddominguez71/deploy-cost-effective-spot-instances-for-testing-environments-in-gcp-with-terraform-7b14a29f2cbf) - IaC approach (for reference, not MVP)

**Confidence Level:** HIGH - Based on official GCP documentation for Spot VMs, startup scripts, and preemption handling. Patterns verified against multiple sources.
