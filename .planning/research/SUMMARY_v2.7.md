# Research Summary: v2.7 Automated GCP Deployment

**Milestone:** v2.7 "Automated GCP Deployment"
**Domain:** Non-interactive cloud deployment automation for scientific computing
**Researched:** 2026-02-06
**Confidence:** HIGH

## Executive Summary

v2.7 replaces interactive v2.6 GCP deployment scripts with a fully automated, config-driven system that finds the cheapest spot instance across all GCP regions and deploys end-to-end without manual intervention. Research reveals this is a well-understood domain with established patterns, but critical pitfalls exist around GCP pricing APIs, Docker resource detection, and spot VM availability.

The recommended approach is **pragmatic automation**: use TOML configuration for human-friendly specs, leverage gcloud CLI (not Terraform or Python APIs), implement a hybrid pricing strategy (CloudPrice.net API with hardcoded regional fallbacks), and dynamically calculate Docker memory/CPU limits on the VM host (not inside containers). This avoids v2.6's critical issues: CPU detection returning 1 inside containers, hardcoded 8GB memory limits on 256GB VMs, and HTTPS failures on bare IP addresses.

Key risks center on **pricing data reliability** (GCP's official Billing API is impractical for spot pricing, third-party APIs may return stale data) and **spot capacity exhaustion** (cheapest region may lack physical hardware despite quota availability). Mitigation: cache pricing with 24-hour TTL, validate outliers, implement multi-region fallback with availability checks, and distinguish quota errors from capacity errors with intelligent retries.

## Key Findings

### Recommended Stack

**Core decision: Extend existing bash infrastructure, don't rewrite.** v2.6's bash scripts are proven and maintainable. The automation layer should augment them with TOML config parsing and pricing queries, not replace them wholesale.

**Core technologies:**
- **TOML config file** (Python 3.11+ tomllib) — Human-readable, stdlib support, Pydantic validation. RESOLVES CONFIG FORMAT DEBATE: TOML wins over YAML (safer, explicit typing, Python ecosystem standard per PEP 621)
- **gcloud CLI non-interactive mode** (--quiet, --format=json) — Official, stable, already in use, supports full automation. Preferred over google-cloud-compute Python library
- **Hybrid pricing strategy** — CloudPrice.net API (daily updates, free) → hardcoded regional rankings (europe-north2, us-central1 cheapest) → fail gracefully. RESOLVES PRICING STRATEGY: Don't use GCP Billing Catalog API directly (too complex, no spot pricing)
- **Pydantic v2** (already in project) — Config validation with clear error messages
- **httpx** (already in project) — CloudPrice API queries with timeout

**Critical version requirements:**
- Python 3.11+ for tomllib stdlib (read-only TOML parser)
- gcloud SDK latest stable (scripting mode support)

### Expected Features

**Must have (table stakes) — from FEATURES.md:**
- Config file with CPU/RAM/disk specs (not machine types)
- Auto-find cheapest spot instance across all regions
- Machine type auto-selection from specs
- Non-interactive execution (zero prompts, works in CI)
- Cost estimate before deployment
- Single-command deployment
- Static IP and persistent disk (already in v2.6)

**Should have (differentiators) — High ROI for v2.7:**
- Multi-region price comparison (not just default region)
- Dry-run mode (--dry-run shows selection without deploying)
- Config validation with pre-flight checks
- Deployment progress logs (timestamped feedback)
- Auto-cleanup on failed deployment (delete orphaned VMs)

**Explicitly NOT building (anti-features):**
- Auto-scaling (single-instance pattern)
- Auto-recovery on preemption (user controls lifecycle)
- Health checks (manual job submission/checking)
- HTTPS/TLS (HTTP-only for fire-and-burn pattern)
- Terraform (overkill for single VM)

**Milestone scope addition:**
- Fix conformer progress tracking display bug (UI issue unrelated to deployment automation)

### Architecture Approach

**Modular bash libraries with single entry point.** Extract reusable functions from v2.6's setup-infrastructure.sh and deploy-vm.sh into lib/ directory. New deploy-auto.sh orchestrates: config load → pricing query → machine selection → infrastructure creation → VM deployment.

**Major components:**
1. **lib/config.sh** — TOML parsing with tomllib, Pydantic validation (CPU/RAM constraints, GCP project existence)
2. **lib/pricing.sh** — CloudPrice.net query with caching, hardcoded fallback, outlier detection
3. **lib/machine.sh** — Map CPU/RAM to machine types via gcloud filtering, calculate Docker memory limits (VM_RAM - 8GB for OS)
4. **lib/infra.sh** — Idempotent infrastructure operations (create-if-missing, reuse-if-exists) extracted from v2.6
5. **deploy-auto.sh** — Orchestrator that sources libraries, handles errors, provides progressive feedback
6. **startup-template.sh** — Dynamic generation based on selected machine specs (NWCHEM_NPROC, WORKER_MEMORY_LIMIT calculated on host)

**Data flow:**
1. config.toml → specs (64 cores, 256GB RAM)
2. Specs → machine types (`gcloud compute machine-types list --filter`)
3. Machine types → pricing (CloudPrice API or hardcoded rankings)
4. Pricing → region selection (cheapest with availability)
5. Specs → Docker limits (WORKER_MEMORY_LIMIT = 248g, NWCHEM_NPROC = 64, written to .env on host)
6. All parameters → deployment (non-interactive gcloud create with --quiet)

**HTTP-only architecture (no HTTPS):**
- No domain validation (IP address only)
- Remove Caddy HTTPS configuration
- Update docker-compose.gcp.yml to expose port 80
- Document security tradeoffs in README

### Critical Pitfalls

**From PITFALLS.md — Top 5 for v2.7:**

1. **GCP Pricing API Complexity (CRITICAL)** — Cloud Billing Catalog API doesn't provide spot pricing directly, requires filtering 500+ SKUs by string matching, and historical data shows stale/incorrect prices (e.g., K80 GPU off by 5x). **Mitigation:** Use CloudPrice.net third-party API with 24-hour cache, validate outliers (if region 90% cheaper than others, it's wrong), fallback to hardcoded regional rankings.

2. **CPU Detection Inside Docker Containers (CRITICAL)** — `nproc` returns 1 inside containers (cgroups isolation), Python's `os.cpu_count()` returns host count but ignores container limits. **v2.6 recurrence:** OpenMPI "not enough slots" error when NWCHEM_NPROC set incorrectly. **Mitigation:** Detect CPUs on VM host BEFORE containerization, write NWCHEM_NPROC to .env file, use --oversubscribe flag (already fixed in v2.6).

3. **Docker Memory Limits Not Matching VM Capacity (CRITICAL)** — v2.6 had 4GB limit on 256GB VM, causing OOM kills mid-calculation. **Mitigation:** Dynamically calculate WORKER_MEMORY_LIMIT = (VM_RAM_GB - 8) on host, write to docker-compose.gcp.yml template. Never use static memory limits.

4. **Zone Resource Exhaustion (MODERATE)** — Cheapest region may lack spot capacity despite quota availability (ZONE_RESOURCE_POOL_EXHAUSTED error). **Mitigation:** Check machine type availability in zone before selection, implement top-3 region fallback list, retry with exponential backoff, distinguish quota errors from capacity errors.

5. **Non-Interactive gcloud Prompts Blocking (CRITICAL)** — gcloud commands wait for user input without --quiet flag, hanging CI/CD pipelines. **Mitigation:** Always use --quiet, set `gcloud config set disable_prompts true`, provide all required parameters explicitly, test in non-TTY environment.

**Known v2.6 issues — prevention checklist:**
- Docker nproc → Fixed, automation must SET NWCHEM_NPROC from VM specs
- Worker memory limit → Fixed with 120GB manual bump, automation must CALCULATE dynamically
- HTTPS on bare IP → Fixed with HTTP, automation must NOT attempt HTTPS
- Conformer progress bug → Unfixed, included in v2.7 scope (separate from deployment)

## Implications for Roadmap

Based on research, suggested phase structure for v2.7:

### Phase 1: Config Foundation and Pricing
**Rationale:** Must establish non-interactive patterns and reliable pricing data before any VM operations. Pricing complexity is CRITICAL pitfall requiring careful validation.

**Delivers:**
- TOML config schema with Pydantic validation
- Config parser (tomllib + validation rules)
- CloudPrice.net pricing query with caching
- Hardcoded regional fallback rankings
- gcloud non-interactive patterns (--quiet, --format=json)
- Unit tests for validation logic

**Addresses features:**
- Config file for compute requirements (table stakes)
- Multi-region price comparison (differentiator)
- Config validation (differentiator)

**Avoids pitfalls:**
- Pitfall #1: GCP Pricing API complexity (implements hybrid strategy)
- Pitfall #5: Non-interactive prompts (establishes --quiet patterns)
- Pitfall #7: API rate limits (caching with 24-hour TTL)

**Research flag:** MEDIUM confidence on CloudPrice API reliability. If API proves unstable during implementation, pivot to hardcoded-only strategy.

### Phase 2: Machine Selection and Resource Calculation
**Rationale:** Must correctly map CPU/RAM specs to machine types and calculate Docker limits BEFORE deployment. CPU/memory detection is CRITICAL pitfall from v2.6.

**Delivers:**
- Machine type filtering via gcloud (CPU/RAM → n2-standard-64 mapping)
- Zone availability checking (prevent ZONE_RESOURCE_POOL_EXHAUSTED)
- Multi-region fallback logic (top 3 cheapest with capacity)
- Dynamic Docker limit calculation (WORKER_MEMORY_LIMIT, NWCHEM_NPROC)
- Startup script template generation with computed values

**Addresses features:**
- Machine type auto-selection (table stakes)
- Auto-find cheapest instance (table stakes)

**Avoids pitfalls:**
- Pitfall #2: CPU detection (calculates on host, not in container)
- Pitfall #3: Memory limits (dynamic calculation, never static)
- Pitfall #4: Zone exhaustion (availability checks + fallback)
- Pitfall #10: Machine type validation (checks availability before selection)

**Research flag:** LOW confidence needed — patterns are well-established, gcloud filtering is documented.

### Phase 3: Deployment Orchestration
**Rationale:** Tie everything together with idempotent infrastructure operations and error handling. Focus on reliability and user experience.

**Delivers:**
- deploy-auto.sh orchestrator sourcing all libraries
- Idempotent infrastructure operations (extracted from v2.6)
- Progressive feedback with timestamps
- Error handling and cleanup on failure
- Cost estimate display before deployment
- Dry-run mode (--dry-run flag)

**Addresses features:**
- Single-command deployment (table stakes)
- Non-interactive execution (table stakes)
- Cost estimate (table stakes)
- Dry-run mode (differentiator)
- Auto-cleanup on failed deploy (differentiator)
- Deployment logs (differentiator)

**Avoids pitfalls:**
- Pitfall #9: Startup script failures (logging to serial console)
- Pitfall #8: Quota vs availability confusion (intelligent error messages)

**Research flag:** LOW confidence needed — orchestration patterns are standard.

### Phase 4: HTTP-Only Container Deployment
**Rationale:** Deploy containers with correct configuration, remove v2.6 HTTPS/Caddy configuration, validate end-to-end.

**Delivers:**
- Updated docker-compose.gcp.yml (HTTP port 80, no Caddy HTTPS)
- Dynamic .env generation with computed limits
- Container startup validation
- Status command updates
- Documentation of HTTP-only pattern

**Addresses features:**
- HTTP-only deployment (anti-feature: NO HTTPS)
- Static IP and persistent disk (already in v2.6)

**Avoids pitfalls:**
- Pitfall #12: HTTPS on bare IP (uses HTTP-only)

**Research flag:** NONE — HTTP deployment is simpler than v2.6 HTTPS attempt.

### Phase 5: Conformer Progress Bug Fix
**Rationale:** Separate UI issue unrelated to deployment automation, can be developed in parallel or after Phase 4.

**Delivers:**
- Fix for conformer progress tracking display
- Updated UI feedback during calculations

**Addresses milestone requirement:**
- Conformer progress bug fix (included in v2.7 scope)

**Research flag:** NONE — UI bug, not related to deployment research.

### Phase Ordering Rationale

**Sequential dependencies:**
- Phase 1 → Phase 2: Pricing data needed before machine selection
- Phase 2 → Phase 3: Machine selection needed before orchestration
- Phase 3 → Phase 4: VM deployment needed before container configuration
- Phase 5 can run in parallel (independent UI fix)

**Risk mitigation order:**
- Phase 1 addresses CRITICAL pricing complexity early (can pivot to hardcoded if API fails)
- Phase 2 prevents v2.6 CPU/memory issues (solved before containers deployed)
- Phase 3 adds reliability patterns on proven foundation
- Phase 4 is simplification (removing HTTPS, lower risk than v2.6)

**Grouping by component:**
- Phase 1+2: Data layer (config, pricing, machine selection)
- Phase 3: Orchestration layer (deployment flow)
- Phase 4: Container layer (Docker configuration)
- Phase 5: UI layer (separate concern)

### Research Flags

**Phases needing deeper research during planning:**
- **Phase 1:** Pricing query implementation — CloudPrice API may require HTML scraping if no JSON endpoint exists, MEDIUM confidence on reliability
- **Phase 2:** Machine type availability checking — Need to validate gcloud machine-types list filtering, but well-documented (LOW risk)

**Phases with standard patterns (skip research-phase):**
- **Phase 3:** Deployment orchestration — Bash scripting patterns, idempotent operations are well-established
- **Phase 4:** Container deployment — Simplification of v2.6 (removing HTTPS), no new patterns
- **Phase 5:** UI bug fix — Not deployment-related, standard frontend debugging

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | TOML resolved over YAML, gcloud CLI over Terraform, Pydantic already in use |
| Features | HIGH | v2.6 experience validates expectations, anti-features clearly defined |
| Architecture | HIGH | Modular bash extraction from v2.6 is low-risk, patterns proven |
| Pitfalls | HIGH | v2.6 post-mortem provides direct experience, GCP docs authoritative |
| Pricing Strategy | MEDIUM | Third-party API dependency is main uncertainty, hardcoded fallback mitigates |

**Overall confidence:** HIGH

v2.7 is an **augmentation**, not a rewrite. Existing v2.6 infrastructure (static IP, persistent disk, firewall, lifecycle scripts) remains unchanged. New automation layer adds config parsing, pricing queries, and dynamic resource calculation — all well-understood patterns. The main uncertainty is CloudPrice.net API reliability, mitigated by hardcoded regional rankings fallback.

### Gaps to Address

**Pricing API reliability (MEDIUM concern):**
- CloudPrice.net may change data format or become unavailable
- **Handling:** Implement early in Phase 1, add fallback to hardcoded rankings, consider GCP Billing Catalog API as tertiary fallback (complex but official)
- **Validation:** Test with stale cache, offline mode, API downtime scenarios

**Machine type availability in practice (LOW concern):**
- Research shows ZONE_RESOURCE_POOL_EXHAUSTED is common, but unclear how OFTEN for 64-core instances
- **Handling:** Implement top-3 region fallback in Phase 2, log which regions failed and why
- **Validation:** Test deployment to exhausted zones, verify fallback behavior

**Conformer progress bug root cause (LOW concern):**
- Separate UI issue, unrelated to deployment automation
- **Handling:** Phase 5 can investigate during implementation
- **Validation:** Not deployment-related, standard debugging

## Milestone Context

**v2.7 replaces v2.6 interactive scripts entirely:**
- v2.6: User prompted for region, zone, machine type, domain
- v2.7: All parameters from config.toml, automated region/type selection, HTTP-only (no domain)

**Key improvements over v2.6:**
- Non-interactive (works in CI/CD, no manual steps)
- Cost-optimized (finds cheapest region automatically)
- Correct resource limits (dynamic calculation, not hardcoded 4GB)
- HTTP-only (removes Caddy HTTPS complexity)
- Dry-run mode (preview before deployment)
- Config validation (fail fast with clear errors)

**Fire-up-and-burn usage pattern:**
- Deploy VM for hours/days, tear down when done
- Manual lifecycle control (no auto-recovery)
- Cost-first optimization (cheapest spot wins)
- HTTP-only (IP address, no domain/HTTPS)
- Ephemeral by design (preemption is acceptable)

## Sources

### Primary (HIGH confidence)

**Official GCP Documentation:**
- [Scripting gcloud CLI](https://cloud.google.com/sdk/docs/scripting-gcloud)
- [Cloud Billing Catalog API](https://docs.cloud.google.com/billing/v1/how-tos/catalog-api)
- [Get Google Cloud pricing information API](https://docs.cloud.google.com/billing/docs/how-to/get-pricing-information-api)
- [Spot VMs | Compute Engine Documentation](https://docs.cloud.google.com/compute/docs/instances/spot)
- [Troubleshooting resource availability](https://docs.cloud.google.com/compute/docs/troubleshooting/troubleshooting-resource-availability)
- [Use startup scripts on Linux VMs](https://cloud.google.com/compute/docs/instances/startup-scripts/linux)
- [gcloud compute machine-types list](https://cloud.google.com/sdk/gcloud/reference/compute/machine-types/list)

**Official Docker Documentation:**
- [Docker Resource Constraints](https://docs.docker.com/engine/containers/resource_constraints/)
- [Docker Compose Deploy Specification](https://docs.docker.com/reference/compose-file/deploy/)

**Python Ecosystem:**
- [Pydantic Settings Management](https://docs.pydantic.dev/latest/concepts/pydantic_settings/)
- [Python tomllib documentation](https://docs.python.org/3/library/tomllib.html)

### Secondary (MEDIUM confidence)

**Third-Party Pricing Tools:**
- [GCP Compute Engine Pricing - CloudPrice](https://cloudprice.net/gcp/compute)
- [GCP Region Pricing Comparison - CloudPrice](https://cloudprice.net/gcp/regions)
- [GCP Preemptible Price History - CloudPrice](https://cloudprice.net/gcp/spot-history)

**Best Practices:**
- [GCP Spot VMs Explained: A Smarter Way to Cut Cloud Costs](https://www.pump.co/blog/spot-instances-gcp)
- [Top 10 GCP cost optimization tools and strategies in 2026](https://northflank.com/blog/gcp-cost-optimization)
- [How to Deploy to GCP in 2026: A Complete Guide](https://encore.cloud/resources/how-to-deploy-to-gcp-2026)
- [JSON vs YAML vs TOML: Which Configuration Format Should You Use in 2026?](https://dev.to/jsontoall_tools/json-vs-yaml-vs-toml-which-configuration-format-should-you-use-in-2026-1hlb)

**Technical Analysis:**
- [Python cpu_count in containers](https://github.com/python/cpython/issues/80235)
- [Docker container nproc issues](https://github.com/moby/moby/issues/43587)
- [CPU count in cgroups](https://donghao.org/2022/01/20/how-to-get-the-number-of-cpu-cores-inside-a-container/)

### Project-Specific (HIGH confidence)

- `.planning/post-v2.6-problems.md` — Direct experience with v2.6 deployment issues (CPU detection, memory limits, HTTPS on bare IP)
- `gcp/` directory — Existing bash scripts to be augmented, not replaced

---

**Research completed:** 2026-02-06
**Ready for requirements:** Yes
**Ready for roadmap:** Yes

**Critical decision resolved:** TOML config (not YAML), hybrid pricing strategy (CloudPrice + hardcoded fallback, not GCP Billing API directly), extend v2.6 bash (not rewrite in Python/Terraform).
