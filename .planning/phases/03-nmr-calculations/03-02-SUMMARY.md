---
phase: 03-nmr-calculations
plan: 02
subsystem: models
tags: [pydantic, api, nmr, models, schemas]

# Dependency graph
requires:
  - phase: 03-nmr-calculations/01
    provides: Preset configurations, solvent validation
provides:
  - AtomShift and NMRResults models for NMR data storage
  - Extended JobInput with preset/solvent fields
  - Extended JobStatusResponse with NMR results
affects: [03-03, calculation pipeline, job storage]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Separate internal models (with shielding) from API response models (shift only)
    - Required solvent field with no default (user must explicitly choose)

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/models.py
    - src/qm_nmr_calc/api/schemas.py
    - src/qm_nmr_calc/api/routers/jobs.py
    - src/qm_nmr_calc/storage.py
    - tests/test_api.py

key-decisions:
  - "Solvent is required field (no default) - users must explicitly choose"
  - "AtomShift stores both shielding and shift; API returns only shift"
  - "NMRResults includes calculation metadata (functional, basis_set, solvent)"

patterns-established:
  - "Internal vs API models: internal has full data, API has user-relevant subset"
  - "Response conversion function handles None checking for optional nested objects"

# Metrics
duration: 4min
completed: 2026-01-19
---

# Phase 3 Plan 2: Data Models for NMR Summary

**Extended Pydantic models and API schemas with preset/solvent input fields and NMR results output**

## Performance

- **Duration:** 4 min
- **Started:** 2026-01-19T18:56:46Z
- **Completed:** 2026-01-19T19:00:16Z
- **Tasks:** 3
- **Files modified:** 5

## Accomplishments
- Added AtomShift and NMRResults models to store NMR calculation results
- Extended JobInput with preset (default: production) and solvent (required) fields
- Extended JobStatus with nmr_results and optimized_geometry_file fields
- Created API response models (AtomShiftResponse, NMRResultsResponse) without internal shielding data
- Extended JobStatusResponse with preset, solvent, and nmr_results fields
- Updated job_status_to_response() to convert NMR results for API output
- Updated storage and router to handle new required solvent field

## Task Commits

Each task was committed atomically:

1. **Task 1: Add NMR result models to models.py** - `8f398ee` (feat)
2. **Task 2: Extend API schemas for NMR** - `4064c72` (feat)
3. **Task 3: Update jobs router conversion function** - `afa55a8` (feat)

## Files Created/Modified
- `src/qm_nmr_calc/models.py` - AtomShift, NMRResults, extended JobInput and JobStatus
- `src/qm_nmr_calc/api/schemas.py` - AtomShiftResponse, NMRResultsResponse, extended JobSubmitRequest and JobStatusResponse
- `src/qm_nmr_calc/api/routers/jobs.py` - job_status_to_response with NMR conversion, endpoints passing new fields
- `src/qm_nmr_calc/storage.py` - create_job_directory accepts solvent and preset parameters
- `tests/test_api.py` - Tests updated to include required solvent parameter

## Decisions Made
- **Solvent required, no default:** Users must explicitly specify the NMR solvent to avoid calculation errors
- **Shielding in internal model only:** Raw isotropic shielding stored for debugging/analysis but not exposed in API
- **NMR results include metadata:** The response includes functional, basis_set, and solvent so users know exactly what calculation was performed

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Storage and tests required updates for new required solvent field**

- **Found during:** Task 3
- **Issue:** The plan stated tests should still pass with new fields, but solvent was marked as required (no default), which broke existing tests and create_job_directory calls
- **Fix:** Updated storage.py to accept solvent parameter, updated jobs router to pass it, updated all tests to include solvent
- **Files modified:** src/qm_nmr_calc/storage.py, src/qm_nmr_calc/api/routers/jobs.py, tests/test_api.py
- **Commit:** afa55a8

## Issues Encountered

None beyond the blocking issue fixed above.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Models ready to store NMR results from calculation pipeline
- API schemas ready to return NMR results to clients
- Storage layer ready to persist jobs with preset/solvent input
- Next plan (03-03) can implement actual NMR calculation and populate these fields

---
*Phase: 03-nmr-calculations*
*Completed: 2026-01-19*
