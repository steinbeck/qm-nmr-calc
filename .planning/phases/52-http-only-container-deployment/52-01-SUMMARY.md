---
phase: 52-http-only-container-deployment
plan: 01
subsystem: infra
tags: [gcp, lifecycle-scripts, toml-config, bash, zone-detection]

# Dependency graph
requires:
  - phase: 51-deployment-orchestration
    provides: TOML config loading via lib/config.sh and deploy-auto.sh orchestrator
  - phase: 50-machine-selection
    provides: Dynamic zone selection based on pricing
provides:
  - All 6 lifecycle scripts (start, stop, delete, status, ssh, logs) using v2.7 TOML config
  - Runtime zone detection via gcloud compute instances list
  - HTTP-only access (no HTTPS/domain/Caddy references)
  - Consistent config loading pattern across all lifecycle scripts
affects: [53-benchmarking-dashboard, deployment-workflows]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Runtime zone detection: Query VM location dynamically instead of hardcoding"
    - "TOML config via lib/config.sh: All scripts source lib/config.sh and call load_config"
    - "HTTP-only deployment: No HTTPS, domain, or Caddy references"

key-files:
  created: []
  modified:
    - gcp/start-vm.sh
    - gcp/stop-vm.sh
    - gcp/ssh-vm.sh
    - gcp/delete-vm.sh
    - gcp/status-vm.sh
    - gcp/logs-vm.sh

key-decisions:
  - "Runtime zone detection: Scripts query gcloud compute instances list to find VM zone dynamically, enabling v2.7's deployment-time zone selection"
  - "HTTP-only URLs: status-vm.sh shows http:// URL instead of domain/HTTPS since v2.7 is HTTP-only"
  - "Config migration: All scripts now load TOML config via lib/config.sh, replacing old config.sh pattern"

patterns-established:
  - "Zone detection pattern: gcloud compute instances list --filter='name=${VM_NAME}' --format='value(zone)'"
  - "Config loading pattern: source lib/config.sh; load_config config.toml"
  - "Lifecycle script independence: Scripts detect VM zone at runtime, work with any deployment location"

# Metrics
duration: 3min
completed: 2026-02-06
---

# Phase 52 Plan 01: HTTP-Only Container Deployment Summary

**All 6 lifecycle scripts migrated to v2.7 TOML config with runtime zone detection, replacing hardcoded zones and HTTPS/domain references with HTTP-only access**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-06T23:43:59Z
- **Completed:** 2026-02-06T23:46:52Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- Migrated all 6 lifecycle scripts from v2.6 config.sh (hardcoded zone) to v2.7 TOML config loading
- Implemented runtime zone detection via gcloud compute instances list in all scripts
- Removed all HTTPS, domain, and Caddy references (v2.7 is HTTP-only)
- Updated all deploy-vm.sh references to deploy-auto.sh
- Scripts now work with VMs deployed to any zone without hardcoding

## Task Commits

Each task was committed atomically:

1. **Task 1: Migrate start-vm.sh, stop-vm.sh, and ssh-vm.sh to v2.7** - `8bc06ce` (feat)
2. **Task 2: Migrate delete-vm.sh, status-vm.sh, and logs-vm.sh to v2.7** - `bc0b390` (feat)

## Files Created/Modified
- `gcp/start-vm.sh` - Load TOML config, detect zone at runtime, show HTTP URL instead of domain
- `gcp/stop-vm.sh` - Load TOML config, detect zone at runtime, reference deploy-auto.sh
- `gcp/ssh-vm.sh` - Load TOML config, detect zone at runtime, reference deploy-auto.sh
- `gcp/delete-vm.sh` - Load TOML config, detect zone at runtime, reference deploy-auto.sh
- `gcp/status-vm.sh` - Load TOML config, detect zone at runtime, remove DOMAIN metadata, show HTTP URL only
- `gcp/logs-vm.sh` - Load TOML config, detect zone at runtime, remove caddy usage reference

## Decisions Made

**1. Runtime zone detection pattern**
- All scripts now query `gcloud compute instances list --filter="name=${VM_NAME}" --format="value(zone)"` to find the VM's zone dynamically
- This enables v2.7's deployment-time zone selection (cheapest zone picked at deploy time)
- Scripts no longer rely on GCP_ZONE from config - zone is discovered at runtime

**2. HTTP-only access URLs**
- status-vm.sh now shows `http://$EXTERNAL_IP` instead of `https://$DOMAIN`
- start-vm.sh updated to show HTTP URL in completion message
- Removed all DOMAIN metadata fetching and display logic

**3. Config loading migration**
- Replaced old config.sh loading block with: `source "$SCRIPT_DIR/lib/config.sh"; load_config "$SCRIPT_DIR/config.toml"`
- All scripts now use TOML config validation via validate_config.py
- Consistent pattern across all lifecycle scripts

**4. deploy-auto.sh references**
- All error messages updated to reference deploy-auto.sh instead of deploy-vm.sh
- Aligns with Phase 51's orchestrator script naming

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all scripts migrated cleanly with consistent patterns.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

All lifecycle scripts ready for Phase 52-02 (teardown-infrastructure.sh migration). Lifecycle scripts now:
- Work with v2.7 TOML config
- Discover VM zone dynamically at runtime
- Show HTTP-only access (no HTTPS/domain)
- Reference deploy-auto.sh consistently

No blockers for next phase.

---
*Phase: 52-http-only-container-deployment*
*Completed: 2026-02-06*
