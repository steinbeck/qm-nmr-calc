---
phase: 27-usage-guide
plan: 02
subsystem: docs
tags: [rest-api, curl, nmr-shifts, documentation, api-reference]

# Dependency graph
requires:
  - phase: 27-01
    provides: Web UI workflow documentation in docs/usage.md
provides:
  - REST API reference with all endpoints documented
  - curl examples for every API endpoint
  - Complete workflow script for programmatic usage
  - Result interpretation guide for shifts, spectra, 3D viewer
  - Ensemble results explanation with Boltzmann averaging
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "RFC 7807 Problem Details for API errors"
    - "Async job submission with 202 Accepted"

key-files:
  created: []
  modified:
    - docs/usage.md

key-decisions:
  - "Documented all endpoints with copy-paste ready curl examples"
  - "Included complete bash script for submit->poll->download workflow"
  - "Explained MAE values and accuracy expectations for different presets"

patterns-established:
  - "API documentation pattern: endpoint, curl example, request/response schemas"
  - "Result interpretation pattern: tables, visualizations, troubleshooting"

# Metrics
duration: 6min
completed: 2026-02-01
---

# Phase 27 Plan 02: REST API Reference and Result Interpretation Summary

**REST API reference with 26 curl examples covering all endpoints, plus comprehensive result interpretation guide for chemical shifts, spectra, and 3D viewer**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-01T07:35:00Z
- **Completed:** 2026-02-01T07:41:00Z
- **Tasks:** 3
- **Files modified:** 1

## Accomplishments

- Documented all REST API endpoints with curl examples (26 examples total)
- Added complete workflow bash script for submit -> poll -> download
- Created result interpretation guide explaining shifts, accuracy, and visualizations
- Documented ensemble results with Boltzmann averaging explanation
- Added troubleshooting section for common result issues

## Task Commits

Each task was committed atomically:

1. **Tasks 1-3: REST API Reference and Result Interpretation** - `e94702d` (docs)

**Note:** Tasks 1-3 were combined into a single commit as they form cohesive documentation additions to the same file.

## Files Created/Modified

- `docs/usage.md` - Added REST API Reference section (~400 lines) and Result Interpretation section (~170 lines)

## Decisions Made

- Used realistic job ID examples (a1b2c3d4e5f6) throughout curl examples
- Included polling pattern with bash while loop for automation
- Provided complete bash script with timeout handling
- Documented all HTTP status codes with examples
- Explained MAE values from DELTA50 benchmark for accuracy expectations

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- docs/usage.md now complete with 944 lines covering:
  - Web UI workflow (from 27-01)
  - Calculation modes and solvent selection
  - REST API reference with all endpoints
  - Result interpretation guide
- Ready for phase 28 (architecture documentation) or 29 (science documentation)

---
*Phase: 27-usage-guide*
*Completed: 2026-02-01*
