---
phase: 50-machine-selection-and-resource-calculation
plan: 02
subsystem: infra
tags: [bash, gcp, machine-selection, wrapper-library, shell-scripting]

# Dependency graph
requires:
  - phase: 50-01
    provides: Python machine selection module with CLI
provides:
  - Bash wrapper library (gcp/lib/machine.sh) for machine selection
  - Four shell functions for VM orchestration scripts
affects: [51-vm-creation, deployment-orchestration, automation-scripts]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Bash wrapper library pattern for Python modules"]

key-files:
  created: [gcp/lib/machine.sh]
  modified: []

key-decisions:
  - "Follow established pattern from pricing.sh and config.sh for consistency"
  - "Use jq for JSON parsing in bash functions"
  - "Provide both JSON output (select_machine) and eval-friendly format (get_docker_resources)"

patterns-established:
  - "Bash wrapper pattern: source library, call functions with args, parse JSON with jq"
  - "script_dir resolution using BASH_SOURCE for reliable path handling"

# Metrics
duration: 43s
completed: 2026-02-06
---

# Phase 50 Plan 02: Machine Selection Bash Wrapper Summary

**Bash wrapper library providing four functions (select_machine, get_docker_resources, generate_startup, get_machine_info) for orchestrating GCP machine selection from shell scripts**

## Performance

- **Duration:** 43 seconds
- **Started:** 2026-02-06T14:47:24Z
- **Completed:** 2026-02-06T14:48:07Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Created gcp/lib/machine.sh with 4 wrapper functions following established pattern
- select_machine() calls Python CLI and returns full JSON output
- get_docker_resources() extracts WORKER_MEMORY_LIMIT and NWCHEM_NPROC for eval pattern
- generate_startup() outputs complete startup script for VM deployment
- get_machine_info() provides human-readable machine type information

## Task Commits

Each task was committed atomically:

1. **Task 1: Create gcp/lib/machine.sh bash wrapper library** - `3a23636` (feat)

**Plan metadata:** (pending - will be committed after SUMMARY and STATE updates)

## Files Created/Modified
- `gcp/lib/machine.sh` - Bash wrapper library with 4 functions wrapping select_machine.py

## Decisions Made

None - followed plan as specified. All wrapper functions implemented following the exact pattern from pricing.sh and config.sh.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None. Implementation was straightforward following established patterns from Phase 49.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Phase 51 (VM Creation and Deployment):**
- Bash wrapper library complete and tested
- All 4 functions verified and working
- Pattern consistent with existing config.sh and pricing.sh wrappers
- Functions provide both JSON output and eval-friendly formats
- 57 tests passing (19 config, 19 pricing, 19 machine)

**What Phase 51 can use:**
- `select_machine 8 32` - get JSON with machine type, zone, region, resources
- `eval $(get_docker_resources 8 32)` - set WORKER_MEMORY_LIMIT and NWCHEM_NPROC vars
- `generate_startup 8 32 qm-nmr-calc 100 > startup.sh` - create startup script
- `get_machine_info 8 32` - display human-readable machine info

**No blockers or concerns.**

---
*Phase: 50-machine-selection-and-resource-calculation*
*Completed: 2026-02-06*
