# Domain Pitfalls: Automated GCP Spot Deployment

**Domain:** Automated GCP spot instance deployment with pricing optimization
**Researched:** 2026-02-06
**Context:** Replacing interactive v2.6 deployment with fully automated, non-interactive system for high-core-count (32-128 cores) NWChem calculations

## Critical Pitfalls

Mistakes that cause deployment failures, data loss, or require complete rewrites.

### Pitfall 1: GCP Pricing API Complexity and Unreliability

**What goes wrong:** The Cloud Billing Catalog API requires separate queries for CPU, RAM, and GPU pricing components, and official pricing data can be dangerously incorrect or lag behind actual prices by days.

**Why it happens:** GCP doesn't use machine families the same way other clouds do. There's no single API endpoint that returns "price for n2-standard-64 spot in us-central1". You must query individual SKUs for CPU type, RAM amount, and combine them. Worse, HTML pricing data at cloud.google.com/compute/all-pricing has been documented as wrong (e.g., Nvidia K80 GPU spot listed as $0.0394/hour but actually $0.19/hour).

**Consequences:**
- Automation selects "cheapest" region based on stale/incorrect data
- Actual costs far exceed estimates
- Pricing email notifications arrive days before OR after price changes, making change detection unreliable

**Prevention:**
1. **Query Cloud Billing Catalog API programmatically**, not HTML scraping
2. **Cache pricing data with TTL** - spot prices change up to once every 30 days, not real-time
3. **Validate pricing results** - if a region is 90% cheaper than others, it's likely wrong
4. **Use pagination correctly** - API returns max 5000 SKUs per page, use nextPageToken
5. **Handle multiprice CUD changes** - as of July 2025, GCP uses more granular per-SKU pricing

**Detection:**
- Price queries return vastly different values for similar regions
- Actual bill doesn't match estimated cost
- API returns incomplete SKU lists (forgot pagination)

**Phase to address:** Phase 1 (pricing query foundation) - build robust SKU resolution and validation before any deployment automation.

**Sources:**
- [GCP Preemptible Price History](https://cloudprice.net/gcp/spot-history)
- [GCP Spot VMs Explained](https://www.pump.co/blog/spot-instances-gcp)
- [Get Google Cloud pricing information API](https://docs.cloud.google.com/billing/docs/how-to/get-pricing-information-api)
- [Cloud Billing Catalog API](https://docs.cloud.google.com/billing/v1/how-tos/catalog-api)
- [GCP CUD changes](https://www.north.cloud/blog/guide-to-google-clouds-multiprice-cud-changes)

### Pitfall 2: CPU Detection Inside Docker Containers

**What goes wrong:** `nproc` returns 1 inside Docker containers, and Python's `os.cpu_count()` returns the HOST CPU count (64 cores) not the container's allocated CPUs (might be 1 or whatever --cpus set).

**Why it happens:** Container isolation via cgroups doesn't update traditional CPU detection tools. `nproc` reads `/proc/cpuinfo` which shows cgroup-isolated view (often 1 core), while `os.cpu_count()` and `multiprocessing.cpu_count()` read host-level information and ignore cgroup limits.

**Consequences:**
- **v2.6 actual issue:** OpenMPI refused to launch 64 processes because `nproc` returned 1, error: "not enough slots available"
- Applications spawn 64 workers thinking they have 64 cores when allocated only 1, causing massive oversubscription
- Memory-per-core calculations wrong (2GB/core × 64 = 128GB when you only have 8GB total)

**Prevention:**
1. **Read cgroup files directly** for accurate CPU limits:
   ```python
   # Read /sys/fs/cgroup/cpu/cpu.cfs_quota_us and cpu.cfs_period_us
   # CPU limit = quota / period (rounded up)
   ```
2. **Use Python's `os.sched_getaffinity(0)`** which respects cgroup affinity (better than `os.cpu_count()`)
3. **For MPI workloads:** Use `--oversubscribe` flag if you KNOW the correct CPU count externally
4. **Set environment variable OUTSIDE container** - detect CPUs on host, pass as env var to container
5. **Never trust `nproc` inside containers**

**Detection:**
- MPI "not enough slots" errors despite requesting high-core VMs
- Container spawns far more workers than allocated CPUs
- CPU utilization metrics show massive oversubscription

**Phase to address:** Phase 2 (VM provisioning) - implement correct CPU detection before any containerized calculations run. This MUST work before moving to Phase 3.

**v2.6 recurrence prevention:** Already fixed with Python `os.cpu_count()` and `--oversubscribe`, but automation must SET NWCHEM_NPROC based on VM specs, not container detection.

**Sources:**
- [Docker container nproc issues](https://github.com/moby/moby/issues/43587)
- [Python cpu_count in containers](https://bugs.python.org/issue36054)
- [os.cpu_count() container limits](https://github.com/python/cpython/issues/80235)
- [CPU count in cgroups](https://donghao.org/2022/01/20/how-to-get-the-number-of-cpu-cores-inside-a-container/)

### Pitfall 3: Docker Memory Limits Not Matching VM Capacity

**What goes wrong:** Docker Compose memory limits set to low defaults (e.g., 4GB) on high-memory VMs (256GB), causing OOM kills or underutilization.

**Why it happens:**
1. **Manual configuration lag** - VM scaled up to 256GB but docker-compose.yml still has `mem_limit: 4g`
2. **Syntax version confusion** - Docker Compose v2 uses `mem_limit`, v3 uses `deploy.resources.limits.memory` (which is IGNORED in non-Swarm mode), leading to limits being silently ignored
3. **Unit parsing errors** - Using "4g" instead of "4294967296" can cause parse failures in some contexts

**Consequences:**
- **v2.6 actual issue:** NWChem with 64 processes needs ~2GB/process = 128GB total, but worker container limited to 4GB caused OOM kills or severe swapping
- Container killed mid-calculation, losing hours of computation
- Automation "succeeds" deploying but calculation fails mysteriously

**Prevention:**
1. **Dynamically set memory limits** based on detected VM RAM:
   ```python
   vm_ram_gb = get_vm_ram()  # Query GCP metadata
   worker_limit = f"{int(vm_ram_gb * 0.9)}g"  # Leave 10% for system
   ```
2. **Use Docker Compose v2.4 syntax** for reliable memory limits (v3 deploy section ignored outside Swarm)
3. **Validate memory limits on startup** - container should check available memory matches expected
4. **Use numeric bytes, not suffixes** if encountering parse errors (4294967296 not "4g")
5. **Test OOM behavior** - ensure container actually gets killed, not just throttled

**Detection:**
- Container killed with exit code 137 (OOM)
- `docker stats` shows memory usage capped well below VM capacity
- Calculation runs but is extremely slow due to swapping

**Phase to address:** Phase 2 (VM provisioning) - template docker-compose.yml with dynamic memory limits based on selected machine type. Do NOT hardcode values.

**v2.6 recurrence prevention:** Automation MUST calculate worker memory limit as function of machine type RAM. Never use static docker-compose files.

**Sources:**
- [Docker memory limits](https://docs.docker.com/engine/containers/resource_constraints/)
- [Docker Compose memory syntax](https://github.com/docker/compose/issues/4513)
- [Memory limit version differences](https://github.com/docker/docs/issues/14185)
- [Docker memory and CPU limits guide](https://oneuptime.com/blog/post/2026-01-16-docker-limit-cpu-memory/view)

### Pitfall 4: Zone/Region Resource Exhaustion

**What goes wrong:** `gcloud compute instances create` fails with ZONE_RESOURCE_POOL_EXHAUSTED or "zone does not have enough resources available" even though you have quota.

**Why it happens:** Spot VMs come from excess Google capacity. High-core-count instances (64-128 cores) are less commonly available. Even if you have quota, the physical hardware might not exist in that zone right now. This is DISTINCT from quota exhaustion - quota is permission, availability is physical capacity.

**Consequences:**
- Deployment script fails after selecting "cheapest" region that has no capacity
- No automatic fallback to second-cheapest region
- User sees error and has to manually retry
- Spot instance creation succeeds but preempted immediately (related but different)

**Prevention:**
1. **Check machine type availability BEFORE selecting region:**
   ```bash
   gcloud compute machine-types list --filter="name=n2-standard-64 AND zone:us-central1-*"
   ```
2. **Implement fallback region list** - top 3-5 cheapest regions, try in order
3. **Use regional MIGs** (Managed Instance Groups) to spread across zones automatically
4. **Select smaller machine types** - GCP explicitly recommends this for spot (64 cores harder than 32 cores)
5. **Retry with exponential backoff** - capacity may become available in minutes
6. **Use `--provisioning-model=SPOT`** explicitly (required for spot instances)

**Detection:**
- Error message contains "ZONE_RESOURCE_POOL_EXHAUSTED" or "not enough resources available"
- Exit code 1 from `gcloud compute instances create`
- Error is NOT quota-related (quota errors are different)

**Phase to address:** Phase 2 (VM provisioning) - implement region selection with availability checking and multi-region fallback. Do NOT assume cheapest region has capacity.

**Sources:**
- [Troubleshooting resource availability](https://docs.cloud.google.com/compute/docs/troubleshooting/troubleshooting-resource-availability)
- [Zone resource pool exhausted](https://bobcares.com/blog/zone_resource_pool_exhausted-gcp-gcp/)
- [Create and use Spot VMs](https://docs.cloud.google.com/compute/docs/instances/create-use-spot)
- [Spot instance availability](https://cast.ai/blog/spot-instance-availability-demystified-aws-azure-and-gcp/)

### Pitfall 5: Non-Interactive gcloud Prompts Blocking Automation

**What goes wrong:** `gcloud` commands wait for user input during automation, causing scripts to hang indefinitely. Prompts like "Do you want to continue (Y/n)?" or "API not enabled, enable it now?" block execution.

**Why it happens:** By default, `gcloud` is interactive. Many commands prompt for confirmation (creating resources, enabling APIs) or ask for missing parameters. In automation/CI, there's no TTY to answer prompts, so command hangs forever.

**Consequences:**
- Deployment script appears stuck, no progress
- CI/CD pipeline times out after 1 hour waiting for prompt
- Silent failures - script waits for input that never comes

**Prevention:**
1. **Always use `--quiet` or `-q` flag** - accepts defaults for all prompts:
   ```bash
   gcloud compute instances create ... --quiet
   ```
2. **Set global config:** `gcloud config set disable_prompts true`
3. **Explicitly provide all required parameters** - don't rely on defaults
4. **Enable APIs beforehand** - don't let gcloud prompt to enable during creation
5. **Check exit codes** - non-zero means error, output may be incomplete
6. **Use `--format` for parseable output** (json/yaml/csv) not human-readable text
7. **Don't depend on stderr messages** - wording changes between gcloud versions

**Detection:**
- Script hangs with no output
- Last log message is a question
- `ps aux` shows `gcloud` process running for hours
- CI job times out at max duration

**Phase to address:** Phase 1 (foundation) - establish gcloud scripting patterns with `--quiet` and full parameter specification. Test in non-interactive environment.

**Sources:**
- [Scripting gcloud CLI](https://cloud.google.com/sdk/docs/scripting-gcloud)
- [Scripting with gcloud guide](https://cloud.google.com/blog/products/management-tools/scripting-with-gcloud-a-beginners-guide-to-automating-gcp-tasks)
- [Non-interactive auth](https://saturncloud.io/blog/authenticate-to-google-container-service-with-script-noninteractive-gcloud-auth-login/)

## Moderate Pitfalls

Mistakes that cause delays, require workarounds, or accumulate technical debt.

### Pitfall 6: Spot VM Preemption During Long Calculations

**What goes wrong:** Spot VM gets preempted by GCP mid-calculation, losing hours of work. NWChem calculations take 30-90 minutes per conformer, but spot instances can be preempted at ANY time with 30 seconds notice.

**Why it happens:** Spot VMs are interruptible by design - GCP reclaims capacity when needed for standard VMs. Preemption probability is "generally low" but varies by zone and day. High-core-count instances may have HIGHER preemption rates (GCP recommends smaller instances for lower rates).

**Consequences:**
- Calculation 80% complete when VM preempted
- No checkpointing = start over from scratch
- Cost savings from spot pricing lost to repeated work

**Prevention:**
1. **Implement checkpointing** - save intermediate results every N minutes
2. **Use preemption handler** - listen for preemption notice, save state within 30 seconds
3. **Consider non-spot for long jobs** - if calculation > 2 hours, evaluate cost vs risk
4. **Regional MIGs with restart** - automatically recreate spot VM after preemption
5. **Break jobs into smaller chunks** - 10× 10-minute tasks better than 1× 100-minute task

**Detection:**
- VM suddenly terminates mid-calculation
- GCP logs show "Instance preempted"
- Calculation restarts from beginning

**Phase to address:** Phase 4 (calculation management) - implement checkpointing and preemption handling. Consider Phase 3 for simple restart automation.

**Sources:**
- [GCP Spot VMs documentation](https://docs.cloud.google.com/compute/docs/instances/spot)
- [Spot VM explained](https://www.pump.co/blog/spot-instances-gcp)

### Pitfall 7: GCP Compute API Rate Limits

**What goes wrong:** Automation hits rate quota limits when querying pricing, checking availability, or creating instances across many regions in quick succession.

**Why it happens:** Compute Engine API has rate quotas (requests per second) that apply per-project. Each `gcloud` command counts toward quota. Aggressive parallel queries (checking all 40+ regions for pricing) can exceed limits.

**Consequences:**
- API returns 429 "Rate limit exceeded"
- Pricing queries fail, automation can't find cheapest region
- Instance creation blocked even though resource quota available

**Prevention:**
1. **Batch API calls** - aggregate requests where possible
2. **Respect rate limits** - space out queries across regions (e.g., 1 request/second)
3. **Cache results** - pricing changes max once per 30 days, cache with TTL
4. **Request quota increase** if legitimate high-volume use
5. **Use exponential backoff** on 429 responses
6. **Prefer aggregated API methods** over individual queries

**Detection:**
- HTTP 429 responses from GCP APIs
- Error: "Rate quota exceeded"
- `gcloud` commands fail with quota errors

**Phase to address:** Phase 1 (pricing queries) - implement rate limiting and caching from the start. Do NOT query all regions in parallel without throttling.

**Sources:**
- [Compute Engine rate quotas](https://cloud.google.com/compute/api-quota)
- [Compute Engine quotas and limits](https://docs.cloud.google.com/compute/quotas-limits)

### Pitfall 8: Quota vs Availability Confusion

**What goes wrong:** Deployment fails with resource unavailability, but team requests quota increase (which doesn't help) instead of trying different zone.

**Why it happens:** Two separate concepts are often confused:
- **Quota:** Permission to use N resources (e.g., 64 CPUs in us-central1)
- **Availability:** Physical capacity exists right now in that zone

You can have quota but no availability. Quota increase doesn't create physical hardware.

**Consequences:**
- Wasted time requesting quota increases that don't solve problem
- User frustration with "quota approved but still failing"

**Prevention:**
1. **Distinguish error types:**
   - Quota: "quota exceeded", exit code 1, mentions specific quota name
   - Availability: "ZONE_RESOURCE_POOL_EXHAUSTED", "not enough resources"
2. **Check quota first:** `gcloud compute project-info describe --project PROJECT_ID`
3. **If quota OK but creation fails:** Availability issue, try different zone
4. **Document difference** in user-facing error messages

**Detection:**
- Resource creation fails but quota checks show plenty available
- Error message says "not enough resources" not "quota exceeded"

**Phase to address:** Phase 2 (VM provisioning) - implement intelligent error handling that distinguishes quota from availability and suggests correct remediation.

**Sources:**
- [Troubleshoot quota errors](https://docs.cloud.google.com/docs/quotas/troubleshoot)
- [Resource availability troubleshooting](https://docs.cloud.google.com/compute/docs/troubleshooting/troubleshooting-resource-availability)

### Pitfall 9: Startup Script Execution Failures

**What goes wrong:** VM created successfully but startup script (that installs Docker, pulls images, starts containers) fails silently. VM exists but application never starts.

**Why it happens:**
- Custom images missing Google Guest Environment (runs startup scripts)
- Script errors with no visible feedback (stderr not logged)
- Script runs as root but needs user context
- Network not ready when script tries to pull Docker images
- Script succeeds but times out (max execution time)

**Consequences:**
- VM appears healthy but ports never open
- Automation waits forever for service to respond
- Silent failure - no clear error message

**Prevention:**
1. **Use public GCP images** (include Guest Environment) or install it manually in custom images
2. **Log everything:** Redirect script output to serial console and Cloud Logging:
   ```bash
   exec > >(tee -a /var/log/startup-script.log)
   exec 2>&1
   ```
3. **Check startup script status:**
   ```bash
   gcloud compute instances get-serial-port-output INSTANCE
   ```
4. **Use metadata key for status** - script writes success/failure to metadata
5. **Add retries for network operations** (apt-get, docker pull)
6. **Test scripts locally** before deploying

**Detection:**
- VM running but application unreachable
- Serial console shows script errors
- Logs show "startup-script exit status 1"

**Phase to address:** Phase 2 (VM provisioning) - build reliable startup script with logging and status reporting before Phase 3 container deployment.

**Sources:**
- [Use startup scripts on Linux VMs](https://cloud.google.com/compute/docs/instances/startup-scripts/linux)
- [Startup script tips and tricks](https://medium.com/google-cloud/few-tips-and-tricks-with-gce-startup-script-323433e2b5ee)
- [Troubleshooting VM startup](https://cloud.google.com/compute/docs/troubleshooting/vm-startup)

### Pitfall 10: Machine Type Selection Without Validation

**What goes wrong:** Automation selects "best" machine type (e.g., c2-standard-60 for 60 cores) but that type doesn't exist in the target region, or spot pricing not available for it.

**Why it happens:** Not all machine types available in all regions. Newer types (C3, N4) limited to specific zones. Some types don't support spot pricing. Automation assumes "if it exists somewhere, it exists everywhere."

**Consequences:**
- Deployment fails after pricing calculation complete
- Error: "machine type not found in zone"
- Have to restart from scratch with different type

**Prevention:**
1. **Check machine type availability in target zone BEFORE selecting:**
   ```bash
   gcloud compute machine-types list \
     --filter="name=c2-standard-60 AND zone:us-central1-a"
   ```
2. **Use machineTypes.aggregatedList API** to map type → available zones
3. **Validate spot support** - not all types support spot pricing
4. **Fallback to similar types** - if c2-standard-60 unavailable, try n2-standard-64
5. **Filter by region early** - don't consider machine types unavailable in target region

**Detection:**
- "Invalid value for field 'machineType'" error
- Machine type name in error message
- Creation fails after region selection

**Phase to address:** Phase 2 (VM provisioning) - integrate machine type availability checking into region/type selection logic. Never assume availability.

**Sources:**
- [Available machine types](https://cloud.google.com/workstations/docs/available-machine-types)
- [Machine families resource guide](https://cloud.google.com/compute/docs/machine-resource)
- [gcloud compute machine-types](https://cloud.google.com/sdk/gcloud/reference/compute/machine-types)
- [View available regions and zones](https://docs.cloud.google.com/compute/docs/regions-zones/viewing-regions-zones)

## Minor Pitfalls

Mistakes that cause annoyance but are easily fixable.

### Pitfall 11: Billing Lag for Cost Tracking

**What goes wrong:** Deployment succeeds, VM runs for 2 hours, but GCP billing shows $0.00. Cost attribution delayed by 1-2 days.

**Why it happens:** GCP doesn't record individual instance costs in real-time. Billing data aggregated and made available to customers 1-2 days later.

**Consequences:**
- Can't provide real-time cost feedback to users
- Budget tracking delayed
- Hard to correlate cost spikes with specific actions

**Prevention:**
1. **Estimate costs locally** based on pricing API + runtime, don't wait for billing
2. **Set billing alerts** for unexpected costs (won't prevent, but alerts after)
3. **Use resource labels** to track which jobs incurred costs (visible when billing data arrives)
4. **Provide estimated cost** to users, not "actual cost from billing"

**Detection:**
- Billing console shows $0 for recently deleted VM
- Cost attribution appears 1-2 days after resource deleted

**Phase to address:** Phase 1 (pricing) or Phase 5 (cost reporting) - implement cost estimation, don't rely on real-time billing data.

**Sources:**
- [GCP Spot VMs Explained](https://www.pump.co/blog/spot-instances-gcp) (mentions billing lag)

### Pitfall 12: HTTPS/TLS on Bare IP Addresses

**What goes wrong:** Automation tries to configure HTTPS/TLS for deployed application, but VM has bare IP (no domain). Can't issue TLS certificate for IP address.

**Why it happens:** Let's Encrypt and similar services don't issue certificates for IP addresses, only domains. Requires domain ownership validation.

**Consequences:**
- **v2.6 actual issue:** Caddy configured for HTTPS but couldn't get cert for 34.44.198.243
- Application unreachable or falls back to HTTP
- User sees certificate errors

**Prevention:**
1. **Use HTTP-only for IP-based deployments** (acceptable for short-lived spot instances)
2. **Use Dynamic DNS** if HTTPS required (map IP to subdomain, get cert for domain)
3. **Use self-signed cert** if HTTPS needed but no domain (user must accept warning)
4. **Cloud Load Balancer** with Google-managed cert (overkill for single VM)

**Detection:**
- Certificate issuance fails
- Let's Encrypt error: "cannot issue for IP addresses"
- Application serves HTTP when HTTPS expected

**Phase to address:** Phase 3 (container deployment) - default to HTTP for spot deployments, document HTTPS limitations.

**v2.6 recurrence prevention:** Already fixed by switching to HTTP. Automation should NOT attempt HTTPS setup.

### Pitfall 13: Platform Differences in Resource Detection

**What goes wrong:** Script works on Linux GCP VM but fails on developer's macOS laptop. Memory detection commands differ (`free` vs `vm_stat`), CPU info locations differ.

**Why it happens:** Development often happens on macOS, deployment on Linux. Commands like `free` (Linux-only) or paths like `/proc/cpuinfo` don't exist on macOS.

**Consequences:**
- **v2.6 actual issue:** Machine info footer feature needed RAM detection, `free` doesn't exist on macOS
- Scripts fail in development, work in production (or vice versa)
- Hard to test locally

**Prevention:**
1. **Use platform-agnostic libraries:** Python's `psutil` works on macOS/Linux/Windows
2. **Platform detection:**
   ```python
   import platform
   if platform.system() == "Darwin":  # macOS
       # use vm_stat
   else:  # Linux
       # use free
   ```
3. **Test in Docker locally** to match production environment
4. **Document platform assumptions** in README

**Detection:**
- Command not found errors on different OS
- Different outputs between macOS and Linux
- Tests pass locally but fail in CI

**Phase to address:** Phase 1 (foundation) - establish platform-agnostic patterns early, especially for resource detection.

**v2.6 recurrence prevention:** Use Python `psutil` or platform detection, never assume Linux-only commands available.

### Pitfall 14: Missing Dependencies in Production

**What goes wrong:** Code works in development (dependency installed locally) but fails in production with ImportError.

**Why it happens:** Dependency used in code but not listed in requirements.txt or Dockerfile. Development environment has it from previous project.

**Consequences:**
- **v2.6 actual issue:** `httpx` used for machine info but not in production dependencies, import failure on deployment
- Application crashes on startup
- "Works on my machine" syndrome

**Prevention:**
1. **Test in clean environment** regularly (fresh Docker build, new virtualenv)
2. **Use dependency scanning:** `pipreqs` to generate requirements.txt from imports
3. **CI runs in isolated environment** catches missing deps before production
4. **Explicit imports check** on app startup, fail fast with clear error

**Detection:**
- ImportError or ModuleNotFoundError in production
- Application works locally but not deployed
- Fresh virtualenv fails

**Phase to address:** Phase 1 (foundation) - establish dependency management practices and CI testing early.

**v2.6 recurrence prevention:** CI should test in clean environment, not rely on cached dependencies.

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| **Phase 1: Pricing Query** | Pricing API complexity, rate limits, pagination | Use Cloud Billing Catalog API with proper pagination, cache results, validate prices |
| **Phase 1: Pricing Query** | Stale/incorrect pricing data | Cache with 24-hour TTL, validate outliers, don't trust HTML pricing |
| **Phase 2: VM Selection** | Machine type unavailable in target region | Check availability before selection, implement fallback types |
| **Phase 2: VM Selection** | Zone resource exhaustion | Multi-region fallback list, check capacity, retry logic |
| **Phase 2: VM Provisioning** | Non-interactive prompts blocking | Always use --quiet, provide all parameters, test in CI |
| **Phase 2: VM Provisioning** | Startup script failures | Log to serial console, check Guest Environment, test scripts |
| **Phase 3: Container Deployment** | Docker memory limits not matching VM | Dynamic memory calculation based on VM specs |
| **Phase 3: Container Deployment** | CPU detection inside containers | Read cgroups or set NWCHEM_NPROC from VM metadata |
| **Phase 3: Container Deployment** | HTTPS on bare IP | Use HTTP-only, document limitation |
| **Phase 4: Calculation Management** | Spot preemption mid-calculation | Checkpointing, preemption handler, or non-spot for long jobs |
| **Phase 4: Cost Tracking** | Billing lag (1-2 days) | Local cost estimation, don't wait for real-time billing |
| **All Phases** | Platform differences (macOS vs Linux) | Use psutil or platform detection, test in Docker |
| **All Phases** | Missing production dependencies | Clean environment testing in CI |

## Known v2.6 Issues - Prevention Checklist

| v2.6 Issue | Status | Automation Must | Phase |
|------------|--------|-----------------|-------|
| Docker nproc returns 1 | Fixed with --oversubscribe | Set NWCHEM_NPROC from VM specs, not container detection | Phase 3 |
| CPU count detection | Fixed with os.cpu_count() | Use Python to detect VM CPUs BEFORE containerization | Phase 2 |
| Worker memory limit too low | Manually bumped to 120GB | Calculate memory limit as f(machine_type.memory_gb) | Phase 2 |
| HTTPS/Caddy on bare IP | Fixed with HTTP | Never attempt HTTPS for spot deployments | Phase 3 |
| Conformer progress tracking bug | Unfixed | Not automation-related, UI issue | N/A |
| macOS vs Linux memory detection | Fixed with platform check | Use psutil or explicit platform branching | Phase 1 |
| Missing httpx dependency | Fixed | CI tests in clean environment | Phase 1 |
| Calculation time expectations | Not a bug | Provide realistic time estimates | Phase 4 |

## Research Confidence

**Overall confidence:** MEDIUM-HIGH

| Area | Level | Reason |
|------|-------|--------|
| GCP Pricing API | MEDIUM | WebSearch + official docs, but API complexity high, conflicting info on reliability |
| CPU/Memory Detection | HIGH | Multiple authoritative sources, GitHub issues with clear explanations, v2.6 real experience |
| Spot Availability | HIGH | Official GCP docs, clear documentation of capacity vs quota |
| gcloud Scripting | HIGH | Official Google scripting guide, well-documented best practices |
| Docker Resource Limits | HIGH | Official Docker docs, GitHub issues with version-specific syntax |
| Startup Scripts | MEDIUM | Official docs but less detail on failure modes |
| v2.6 Issue Prevention | HIGH | Direct project experience, known root causes |

## Sources Summary

**Authoritative (Official GCP/Docker Documentation):**
- [Scripting gcloud CLI](https://cloud.google.com/sdk/docs/scripting-gcloud)
- [Cloud Billing Catalog API](https://docs.cloud.google.com/billing/v1/how-tos/catalog-api)
- [Troubleshooting resource availability](https://docs.cloud.google.com/compute/docs/troubleshooting/troubleshooting-resource-availability)
- [Compute Engine rate quotas](https://cloud.google.com/compute/api-quota)
- [Docker resource constraints](https://docs.docker.com/engine/containers/resource_constraints/)
- [GCP Spot VMs documentation](https://docs.cloud.google.com/compute/docs/instances/spot)
- [Use startup scripts on Linux VMs](https://cloud.google.com/compute/docs/instances/startup-scripts/linux)

**Community/Technical Analysis:**
- [GCP Spot VMs Explained](https://www.pump.co/blog/spot-instances-gcp)
- [Spot instance availability analysis](https://cast.ai/blog/spot-instance-availability-demystified-aws-azure-and-gcp/)
- [Python cpu_count() container issue](https://github.com/python/cpython/issues/80235)
- [Docker nproc issue](https://github.com/moby/moby/issues/43587)
- [Startup script tips](https://medium.com/google-cloud/few-tips-and-tricks-with-gce-startup-script-323433e2b5ee)

**Project-Specific:**
- `.planning/post-v2.6-problems.md` (direct experience with v2.6 deployment)
