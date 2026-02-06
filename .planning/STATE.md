# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-06)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.7 Automated GCP Deployment

## Current Position

Milestone: v2.7 Automated GCP Deployment
Phase: Phase 52 - HTTP-Only Container Deployment (COMPLETE, VERIFIED)
Plan: 2/2 complete, 14/14 must-haves verified
Status: Phase 52 complete and verified, ready for Phase 53
Last activity: 2026-02-06 — Phase 52 verified (14/14 must-haves, 4/4 requirements satisfied)

Progress: [#####################] 112 plans complete (v1.0-v2.6 + 49-01, 49-02, 50-01, 50-02, 51-01, 51-02, 52-01, 52-02)

## Performance Metrics

**Velocity:**
- Total plans completed: 112 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10, v2.3: 3, v2.4: 8, v2.5: 4, v2.6: 4, v2.7: 8)
- Average duration: ~6.4 min
- Total execution time: ~715 min (~11.9 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration | Status |
|-----------|--------|-------|----------|--------|
| v1.0 Core NMR Service | 6 | 16 | 2 days | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | ~2 days | Shipped 2026-01-28 |
| v2.0.1 Conformer Pre-selection | 1 | 3 | 18 min | Shipped 2026-01-30 |
| v2.1 UI Redesign | 6 | 17 | ~3 days | Shipped 2026-01-31 |
| v2.2 Documentation | 7 | 10 | 2 days | Shipped 2026-02-01 |
| v2.3 NMReData Export | 3 | 3 | 1 day | Shipped 2026-02-01 |
| v2.4 Docker Deployment | 6 | 8 | ~2 hours | Shipped 2026-02-03 |
| v2.5 ARM64 Docker Support | 4 | 4 | ~1 day | Shipped 2026-02-04 |
| v2.6 GCP Spot Deployment | 5 | 4 | ~1 day | Shipped 2026-02-05 |
| v2.7 Automated GCP Deployment | 5 | 8/TBD | In progress | Started 2026-02-06 |

## Accumulated Context

### Decisions

All prior decisions logged in PROJECT.md Key Decisions table.
v2.6 decisions archived.

**v2.6 Deployment Problems (see .planning/post-v2.6-problems.md):**
- Docker nproc returns 1 inside containers — fixed with --oversubscribe
- Worker memory limit too low for high-memory VMs — needed manual bumps
- HTTPS/Caddy fails on bare IP — had to switch to HTTP
- Conformer progress tracking display bug — still unfixed, addressed in v2.7 Phase 53
- macOS vs Linux memory detection differences
- Missing httpx production dependency

**v2.7 Approach (from research):**
- TOML config (not YAML) for compute requirements
- CloudPrice.net API with hardcoded fallback (not GCP Billing API directly)
- Extend v2.6 bash scripts with automation libraries (not rewrite in Python/Terraform)
- Dynamic Docker limit calculation on VM host (not inside containers)
- HTTP-only deployment (no domain, no HTTPS)
- Modular bash libraries (config.sh, pricing.sh, machine.sh, infra.sh)

**v2.7 Phase 50 (Machine Selection) decisions:**
- find_available_zone() calls get_ranked_regions() internally for simpler API
- Docker worker memory = VM RAM - 8GB OS overhead (minimum 4GB)
- Startup script uses runtime $(nproc) detection instead of hardcoded CPU count
- HTTP-only docker-compose.gcp.yml override (no Caddy service)
- --oversubscribe flag for MPI compatibility in containers (v2.6 fix)
- Bash wrapper library pattern: source gcp/lib/machine.sh, call functions with args
- Provide both JSON output (select_machine) and eval-friendly format (get_docker_resources)

**v2.7 Phase 51 (Deployment Orchestration) decisions:**
- Cleanup trap only deletes VMs, never disks or IPs (data preservation)
- Trap registration deferred to orchestrator, not set in library (composability)
- Existence checks always run, not wrapped in execute() (accurate dry-run reporting)
- HTTP and SSH firewall rules only, no HTTPS (v2.7 HTTP-only per RPL-03)
- Orchestrator script replaces v2.6 deploy-vm.sh completely
- Startup script generated to temp file and cleaned up after VM creation

**v2.7 Phase 52 (HTTP-Only Container Deployment) decisions:**
- Runtime zone detection: Scripts query gcloud compute instances list to find VM zone dynamically
- HTTP-only URLs: status-vm.sh shows http:// URL instead of domain/HTTPS
- Config migration: All scripts now load TOML config via lib/config.sh
- deploy-auto.sh references: All error messages updated to reference deploy-auto.sh
- Teardown zone detection: try VM first, fallback to disk, handle neither existing gracefully
- Region derivation from zone: ${GCP_ZONE%-*} to strip zone suffix
- Conditional cleanup: disk/IP deletion based on zone availability with fallback strategies

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), shipped 2026-01-28
- v2.0.1: 1 phase (24), shipped 2026-01-30
- v2.1: 6 phases (18-23), shipped 2026-01-31
- v2.2: 7 phases (25-31), shipped 2026-02-01
- v2.3: 3 phases (32-34), shipped 2026-02-01
- v2.4: 6 phases (35-40), shipped 2026-02-03
- v2.5: 4 phases (41-44), shipped 2026-02-04
- v2.6: 5 phases (45-48.1), shipped 2026-02-05
- v2.7: 5 phases (49-53), roadmap complete 2026-02-06

### Pending Todos

- Dark mode (deferred from v2.1)
- Enhanced interactivity (deferred from v2.1)
- User accounts and calculation history
- NMReData enhancements: INCHI tag, 1D pseudo-spectrum tags, J-coupling
- Phase 48.1 implementation (machine info display - v2.6)

### Blockers/Concerns

**Active:**
None

**Resolved (v2.7 planning):**
- Pricing strategy: Use CloudPrice.net API with hardcoded fallback
- Config format: TOML selected over YAML
- Implementation approach: Augment v2.6 bash scripts, don't rewrite

## Session Continuity

Last session: 2026-02-06
Stopped at: Completed 52-02-PLAN.md (teardown infrastructure migration)
Resume file: None
Next: Phase 53 (conformer progress bug fix)
Tests: 415 tests collected (377 pre-existing + 19 config + 19 pricing + 19 machine; pre-existing failures in NWChem integration unrelated)
Codebase: ~7,300 LOC Python, ~3,050 LOC tests, ~950 LOC templates, ~2,400 LOC CSS, ~4,800 LOC docs, ~2,170 LOC GCP scripts
Docker: Worker image 2.1GB (amd64), API image ~733MB (multi-arch), ARM64 worker 2.1GB (arm64), multi-arch manifests on GHCR
GCP: Phase 52 complete - all 7 lifecycle/teardown scripts migrated to v2.7 TOML config + runtime zone detection, HTTP-only deployment verified
