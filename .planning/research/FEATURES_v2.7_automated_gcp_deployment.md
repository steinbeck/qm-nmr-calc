# Feature Landscape: Non-Interactive GCP Spot Deployment

**Domain:** Fire-up-and-burn compute deployment (ephemeral scientific workloads)
**Researched:** 2026-02-06
**Confidence:** HIGH

## Executive Summary

Non-interactive GCP spot deployment tools for scientific computing must balance **cost optimization** (finding cheapest instances) with **zero-friction deployment** (no manual steps). This domain has clear table stakes (config files, region/zone selection, persistent data) and equally clear anti-features (auto-scaling, health checks, auto-recovery).

The fire-up-and-burn pattern is fundamentally different from production services:
- **Ephemeral by design** - Deploy for hours/days, tear down when done
- **Manual lifecycle** - User controls start/stop, no auto-recovery
- **Cost-first** - Cheapest spot instance wins, availability is secondary
- **HTTP-only** - No domains, no HTTPS, just IP addresses

## Table Stakes

Features users expect. Missing = product feels incomplete.

| Feature | Why Expected | Complexity | Implementation Notes |
|---------|--------------|------------|---------------------|
| **Config file for compute requirements** | Can't be interactive if user has to answer prompts | Low | YAML/JSON with `cpu_count`, `memory_gb`, `disk_gb` |
| **Auto-find cheapest spot instance** | Manual price comparison defeats "fire-up-and-burn" value prop | Medium | Query GCP pricing API, filter by region/zone/availability |
| **Machine type auto-selection** | Users specify resources, not machine types | Medium | Map CPU/RAM to machine families (n2-standard, c2-standard, etc.) |
| **Static IP reservation** | Needed for DNS, persistent across VM recreation | Low | Reserve once during setup, reuse forever |
| **Persistent disk for data** | Job data must survive VM deletion | Low | Separate disk lifecycle from VM lifecycle |
| **Graceful shutdown on preemption** | Spot VMs get 30s notice, must handle cleanly | Medium | Shutdown script with Docker Compose down, state persistence |
| **Single-command deployment** | "Deploy and go" - no multi-step wizards | Low | One script orchestrates everything |
| **Cost estimate before deployment** | Users need to approve spend | Low | Calculate hourly cost from machine type + disk + IP |
| **Non-interactive execution** | Must work in CI/CD, scripted environments | Medium | All config from file, no stdin prompts |
| **Lifecycle scripts (start/stop/delete)** | Manual control over when to spend money | Low | Thin wrappers around `gcloud compute instances start/stop/delete` |
| **Status command** | Quick check: "Is it running? What's the IP?" | Low | `gcloud compute instances describe` with formatted output |

**v2.6 already has:** Static IP, persistent disk, lifecycle scripts, status command, cost estimates, graceful shutdown

**Missing for v2.7:** Config file, auto-find cheapest instance, machine type auto-selection, non-interactive execution

## Differentiators

Features that make it really smooth. Not expected, but valued.

| Feature | Value Proposition | Complexity | Implementation Notes |
|---------|-------------------|------------|---------------------|
| **Multi-region price comparison** | Find absolute cheapest across ALL GCP regions | Medium | Parallel query all regions, rank by price, warn about latency |
| **Availability-aware selection** | Avoid regions with frequent preemption | Medium | Query GCP spot availability history, deprioritize low-availability zones |
| **Dry-run mode** | See what would be deployed without deploying | Low | `--dry-run` flag shows selected machine, region, cost |
| **SSH key auto-injection** | No manual key management | Low | Use `gcloud compute config-ssh` or startup script |
| **Docker image pre-pull** | Faster startup (images pulled before "ready" signal) | Low | Startup script pulls GHCR images during provisioning |
| **Config validation** | Catch errors before deployment | Low | Pre-flight checks: valid CPU/RAM combos, disk size, quotas |
| **Deployment logs** | See what's happening during startup | Medium | Stream startup script logs to local terminal |
| **Auto-cleanup on failed deploy** | Don't leave orphaned VMs if deployment fails | Medium | Trap errors, delete VM if startup script fails |
| **Region blacklist/whitelist** | "Never use asia-east, always prefer us-central" | Low | Config file: `allowed_regions`, `blocked_regions` |
| **Spot price alerts** | Notify if spot price spikes | Low | Cache selected price, warn if current price > 2x cached |
| **One-command teardown** | Delete everything (VM + infrastructure) | Low | `teardown-all.sh` removes VM, IP, disk, firewall |
| **Resume from preemption** | Re-deploy to same config after preemption | Medium | Cache last deployment config, `--resume` flag |

**High-value differentiators for v2.7:**
1. **Multi-region price comparison** - Core value prop
2. **Dry-run mode** - Builds confidence before spend
3. **Config validation** - Saves debugging time
4. **Deployment logs** - Visibility during startup
5. **Auto-cleanup on failed deploy** - Prevents orphaned resources

## Anti-Features

Features to explicitly NOT build. Common mistakes in this domain.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Auto-scaling** | Fire-up-and-burn is single-instance, fixed size | User requests N cores, gets N cores, done |
| **Auto-recovery on preemption** | Preemption = expected, manual re-run is fine | Log preemption event, user runs `deploy-vm.sh` again |
| **Health checks** | No need - user submits job, checks results manually | Provide status command, don't poll |
| **Load balancers** | Single VM deployment, no multi-instance traffic | Direct IP access is fine |
| **Managed instance groups (MIGs)** | Adds complexity for zero benefit in single-instance case | Use standalone VM |
| **Kubernetes/GKE** | Massive overkill for "run Docker Compose on one VM" | Docker Compose is table stakes here |
| **HTTPS/TLS certificates** | HTTP-only pattern, no domain names | Serve on port 80, document security tradeoffs |
| **Monitoring/metrics** | Short-lived workloads don't need Prometheus | Status command + container logs is enough |
| **Rolling updates** | No downtime requirement, stop-deploy-start is fine | Full teardown/recreate pattern |
| **Persistent sessions** | Job state in database, not user sessions | Stateless HTTP API |
| **Multi-tenancy** | One user, one VM, disposable | Single deployment per user |

**Critical anti-feature for v2.7:** Do NOT build auto-recovery. Preemption should be logged and visible, but user decides when to re-run.

## Feature Dependencies

```
Non-interactive deployment requires:
├─ Config file (CPU, RAM, disk) ─────────┐
├─ GCP pricing query ────────────────────┤
├─ Machine type selection ───────────────┤──> Cheapest instance selection
├─ Region/zone filtering ────────────────┘
│
├─ Static IP (from v2.6) ────────────────┐
├─ Persistent disk (from v2.6) ──────────┤──> VM deployment
├─ Firewall rules (from v2.6) ───────────┤
├─ Startup script (from v2.6) ───────────┤
└─ Docker Compose (from v2.6) ───────────┘

Dry-run mode requires:
└─ All selection logic (without actual deployment)

Multi-region comparison requires:
├─ GCP pricing API access
└─ Parallel region queries (performance optimization)

Resume from preemption requires:
├─ Last deployment config cache
└─ Idempotent deployment script
```

## MVP Recommendation

For v2.7 MVP, prioritize:

### Must Have (Table Stakes)
1. **Config file** - `gcp/compute-config.yml` with CPU/RAM/disk
2. **Machine type auto-selection** - Map config to machine families
3. **GCP pricing query** - Get spot prices for selected regions
4. **Cheapest instance selection** - Rank by price, pick #1
5. **Non-interactive deployment** - Zero prompts, all from config
6. **Cost estimate display** - Show before deployment, require confirmation via flag

### Should Have (High-Value Differentiators)
7. **Multi-region price comparison** - Check all regions, not just default
8. **Dry-run mode** - `--dry-run` shows selection without deploying
9. **Config validation** - Pre-flight checks before API calls
10. **Deployment logs** - Stream startup progress to terminal

### Nice to Have (Defer to post-v2.7)
- Region blacklist/whitelist
- Availability-aware selection
- Spot price alerts
- Resume from preemption
- Auto-cleanup on failed deploy

## Defer to Post-MVP

**Good ideas, wrong milestone:**
- **Terraform/IaC support** - Shell scripts work fine for MVP
- **Web UI for config** - Config file is simpler
- **Cost budgets/alerts** - User controls lifecycle manually
- **GPU instance support** - Different pricing/availability logic
- **Snapshot/backup automation** - Persistent disk already handles this
- **Custom machine types** - Adds API complexity, standard types cover 90% of use cases

## v2.6 → v2.7 Migration

**What v2.6 already provides (keep):**
- Infrastructure setup scripts (`setup-infrastructure.sh`)
- Lifecycle management (`start-vm.sh`, `stop-vm.sh`, `delete-vm.sh`, `status-vm.sh`)
- Persistent disk for data
- Static IP reservation
- Firewall rules (HTTP, SSH)
- Startup script with Docker installation
- Shutdown script with graceful container stop
- Docker Compose GCP override file
- SSH access and log streaming

**What v2.7 replaces:**
- Interactive prompts in `deploy-vm.sh` → Config file
- Manual region/zone selection → Auto-selection based on price
- Manual machine type selection → Auto-selection based on CPU/RAM requirements
- DNS validation → Remove (HTTP-only pattern)
- Domain prompts → Remove (HTTP-only pattern)

**What v2.7 adds:**
- GCP pricing API integration
- Machine type database/mapping
- Multi-region price comparison
- Dry-run mode
- Config file parser
- Pre-flight validation

## Implementation Patterns from Research

### GCP Pricing Optimization (Sources: [1][2][3])

**Best practice in 2026:** Use `gcloud compute machine-types list` to enumerate options, cross-reference with spot pricing history API.

**Spot pricing volatility:** Spot prices change up to once per day, but are stable enough for "deploy now" decisions. Don't optimize for future price predictions.

**Region selection:** Smaller machine types have better spot availability. Multi-zone distribution increases availability but adds complexity (anti-pattern for single-VM deployment).

**Cost savings:** Spot VMs provide 60-91% discount vs on-demand, making pricing optimization the #1 feature for cost-conscious users.

### Config File Patterns (Sources: [4][5])

**YAML over JSON:** Better for human editing, comments supported. Example:
```yaml
compute:
  cpu_cores: 64
  memory_gb: 128
  disk_gb: 100

deployment:
  allowed_regions:
    - us-central1
    - us-east1
  prefer_region: us-central1  # Latency optimization
```

**Validation before deployment:** Check CPU/RAM combinations are valid (e.g., n2-standard requires 4GB per core).

### Graceful Shutdown (Sources: [6][7])

**30-second window:** GCP Spot VMs provide 30s notice before preemption. For Kubernetes/GKE, grace period is only 15s, but standalone VMs get full 30s.

**Shutdown script best practices:**
1. Trap SIGTERM in startup script
2. Run `docker compose down` (20-25s is enough)
3. Flush logs to persistent disk
4. Exit cleanly

**Current v2.6 implementation:** Already has shutdown.sh with graceful stop. Keep as-is.

### Fire-and-Forget Patterns (Sources: [8][9])

**Ephemeral infrastructure:** Resources created dynamically, destroyed when not needed. This is the defining characteristic of fire-up-and-burn.

**Security considerations:** Ephemeral workloads skip runtime security scanning (too slow for short lifecycles). Apply security at image build time instead.

**Cost optimization:** Automatic shutdown when idle is a key pattern, but requires job completion detection (defer to post-MVP).

### Machine Type Selection (Sources: [10][11])

**Attribute-based selection:** Specify vCPUs, memory, storage as attributes, let GCP find matching machine types. This is the recommended 2026 pattern.

**Implementation:** Use `gcloud compute machine-types list --filter="zone:(us-central1-a)"` to get available types, filter by CPU/RAM requirements, rank by spot price.

**Machine families:**
- `n2-standard` - Balanced CPU/RAM (4GB per core)
- `n2-highmem` - Memory-intensive (8GB per core)
- `n2-highcpu` - CPU-intensive (1GB per core)
- `c2-standard` - Compute-optimized
- `e2-*` - Cheapest, limited performance

**For DFT calculations:** Prefer `n2-standard` or `c2-standard` (compute-bound workload).

## Research Confidence Assessment

| Topic | Confidence | Rationale |
|-------|------------|-----------|
| Table stakes features | HIGH | v2.6 implementation validates user expectations |
| GCP pricing optimization | HIGH | Official documentation + 2026 best practices |
| Config file patterns | HIGH | Industry standard for IaC tools |
| Graceful shutdown | HIGH | Official GCP documentation, tested in v2.6 |
| Fire-and-forget pattern | MEDIUM | Validated by research, but less specific to GCP |
| Machine type selection | HIGH | GCP attribute-based selection is current (2026) best practice |
| Anti-features | HIGH | v2.6 experience confirms what NOT to build |

## Sources

1. [GCP Spot VMs Explained: A Smarter Way to Cut Cloud Costs](https://www.pump.co/blog/spot-instances-gcp)
2. [Top 10 GCP cost optimization tools and strategies in 2026](https://northflank.com/blog/gcp-cost-optimization)
3. [Spot VMs Pricing | Google Cloud](https://cloud.google.com/spot-vms/pricing)
4. [How to Deploy to GCP in 2026: A Complete Guide](https://encore.cloud/resources/how-to-deploy-to-gcp-2026)
5. [Understanding configurations | Cloud Deployment Manager](https://cloud.google.com/deployment-manager/docs/step-by-step-guide/create-a-configuration)
6. [Spot VMs | Compute Engine Documentation](https://docs.cloud.google.com/compute/docs/instances/spot)
7. [Spot VMs | GKE Documentation](https://docs.cloud.google.com/kubernetes-engine/docs/concepts/spot-vms)
8. [11 cloud cost optimization strategies for 2026](https://northflank.com/blog/cloud-cost-optimization)
9. [Ephemeral Workloads Security in Cloud Environments](https://www.gitguardian.com/nhi-hub/ephemeral-workload-security-in-cloud-environments)
10. [Amazon EC2 Spot Instances Guide 2026](https://sedai.io/blog/optimizing-spot-instances-in-aws)
11. [Future Proof Cost Optimization with Attribute-Based Instance Type Selection](https://aws.amazon.com/blogs/apn/future-proof-cost-optimization-with-attribute-based-instance-type-selection-and-amazon-ec2-spot/)
