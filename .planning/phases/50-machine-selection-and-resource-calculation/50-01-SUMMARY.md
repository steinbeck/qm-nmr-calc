---
phase: 50-machine-selection-and-resource-calculation
plan: 01
subsystem: infra
tags: [gcp, compute-engine, machine-types, docker, resource-limits, startup-scripts]

# Dependency graph
requires:
  - phase: 49-config-foundation-and-pricing-query
    provides: GCPConfig model, get_ranked_regions() pricing API with cache/fallback
provides:
  - Machine type selection via gcloud with CPU/RAM filtering
  - Zone validation with regional fallback using get_ranked_regions()
  - Docker resource calculation (VM RAM - 8GB OS overhead)
  - HTTP-only startup script generator with nproc detection
  - CLI for machine selection and startup script generation
affects: [51-vm-creation-and-deployment, 52-automation-orchestration]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - TDD with 19 test cases all passing (RED-GREEN cycle)
    - Single gcloud wrapper (_run_gcloud) as mock target for all tests
    - Integration with get_ranked_regions() for zone selection
    - HTTP-only deployment (no Caddy/HTTPS)

key-files:
  created:
    - gcp/select_machine.py
    - tests/test_gcp_machine.py
  modified: []

key-decisions:
  - "find_available_zone() calls get_ranked_regions(cpu_cores, ram_gb) internally rather than receiving regions as parameter"
  - "Docker worker memory limit = VM RAM - 8GB OS overhead (minimum 4GB after overhead)"
  - "Startup script uses runtime nproc detection instead of hardcoded CPU count"
  - "HTTP-only deployment on port 80 (no Caddy service)"
  - "--oversubscribe flag in MPI command for Docker container compatibility (v2.6 fix)"

patterns-established:
  - "TDD pattern: Create all tests first (RED), implement until all pass (GREEN), commit atomically"
  - "Single mock target pattern: _run_gcloud() wrapper isolates all gcloud calls for easy testing"
  - "Integration pattern: find_available_zone() transparently calls get_ranked_regions() from query_pricing.py"

# Metrics
duration: 4min
completed: 2026-02-06
---

# Phase 50 Plan 01: Machine Selection and Resource Calculation Summary

**Machine type selection with zone fallback via get_ranked_regions(), Docker resource limits (VM RAM - 8GB), and HTTP-only startup script generation with nproc detection**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-06T15:14:10Z
- **Completed:** 2026-02-06T15:18:23Z
- **Tasks:** 2 (TDD: RED + GREEN)
- **Files modified:** 2

## Accomplishments

- Machine type selection maps CPU/RAM requirements to smallest matching GCP machine via gcloud
- Zone validation with regional fallback using get_ranked_regions() integration
- Docker resource calculation with 8GB OS overhead, validates minimum 4GB available
- HTTP-only startup script generator with runtime nproc detection
- Complete CLI with JSON output and --generate-startup-script modes
- All 19 tests pass, no regressions (57/57 GCP tests pass)

## Task Commits

Each TDD phase was committed atomically:

1. **Task 1: Write failing tests (RED)** - `b3ea0da` (test)
   - 19 test cases covering all 5 functions + CLI
   - All fail with ImportError (module doesn't exist yet)

2. **Task 2: Implement module (GREEN)** - `8136df9` (feat)
   - select_machine_type, validate_machine_in_zone, find_available_zone
   - calculate_docker_resources, generate_startup_script
   - CLI with argparse supporting both JSON and script output
   - All 19 tests pass

**Plan metadata:** Not yet committed (pending in this summary update)

## Files Created/Modified

- `gcp/select_machine.py` (594 lines) - Machine selection, validation, resource calculation, startup script generation with CLI
- `tests/test_gcp_machine.py` (253 lines) - 19 test cases with mocked gcloud calls

## Decisions Made

**1. Internal get_ranked_regions() integration**
- find_available_zone() calls get_ranked_regions(cpu_cores, ram_gb) internally
- Caller doesn't need to know about pricing module
- Simpler API: just pass CPU/RAM requirements, get back best zone

**2. Docker memory overhead calculation**
- Worker memory limit = VM RAM - 8GB OS overhead
- Validates minimum 4GB available after overhead
- Format: "{N}g" suffix for Docker memory limits

**3. Runtime nproc detection**
- Startup script uses $(nproc) at runtime instead of hardcoded CPU count
- Handles VM types where gcloud CPU count might differ from actual
- More robust for different GCP machine families

**4. HTTP-only deployment**
- No Caddy service in docker-compose.gcp.yml
- Direct port 80 exposure on API service
- Simplified for spot instances without domains

**5. MPI --oversubscribe flag**
- Carried forward from v2.6 fix
- Required for MPI to work inside Docker containers
- Avoids "not enough slots available" errors

## Deviations from Plan

**1. [Rule 1 - Bug] Fixed bash variable escaping in startup script**
- **Found during:** Task 2 (GREEN phase)
- **Issue:** Python f-string had `\$` which caused SyntaxWarning (invalid escape sequence)
- **Fix:** Removed backslash escapes for bash variables inside heredoc (bash handles escaping)
- **Files modified:** gcp/select_machine.py
- **Verification:** All tests pass with no warnings
- **Committed in:** 8136df9 (part of Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Minor syntax fix for Python warnings. No functional change.

## Issues Encountered

None - TDD execution was smooth. All tests passed on first GREEN attempt after syntax fix.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

Ready for Phase 51 (VM Creation and Deployment):
- Machine type selection complete
- Zone fallback logic tested
- Docker resource limits calculated
- Startup script generator ready

All functions export properly and integrate with existing GCP modules (validate_config.py, query_pricing.py).

---
*Phase: 50-machine-selection-and-resource-calculation*
*Completed: 2026-02-06*
