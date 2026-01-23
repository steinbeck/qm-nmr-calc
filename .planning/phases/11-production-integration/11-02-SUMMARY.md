---
phase: 11-production-integration
plan: 02
subsystem: models
tags: [pydantic, backwards-compatibility, storage, api]

# Dependency graph
requires:
  - phase: 11-production-integration
    plan: 01
    provides: "ISiCLE dependency removed from project"
provides:
  - "JobStatus model without isicle_version field"
  - "Backwards compatibility for old job files with isicle_version"
  - "Storage and API simplified without isicle_version parameters"
affects: [job-creation, api-responses, storage]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Pydantic extra='ignore' for backwards-compatible field removal"
    - "Model evolution strategy for removing deprecated fields"

key-files:
  created: []
  modified:
    - "src/qm_nmr_calc/models.py"
    - "src/qm_nmr_calc/storage.py"
    - "src/qm_nmr_calc/api/routers/jobs.py"
    - "src/qm_nmr_calc/api/routers/web.py"

key-decisions:
  - "Use extra='ignore' in model_config to allow old job files with isicle_version to load without error"
  - "Remove isicle_version parameter completely from create_job_directory() signature"

patterns-established:
  - "Backwards-compatible field removal: extra='ignore' + remove field definition"
  - "Simplified function signatures after deprecated parameters removed"

# Metrics
duration: 6min
completed: 2026-01-23
---

# Phase 11 Plan 02: ISiCLE Version Removal Summary

**JobStatus model and storage layer cleaned of ISiCLE references while maintaining backwards compatibility for old job files via Pydantic extra='ignore'**

## Performance

- **Duration:** 6 min
- **Started:** 2026-01-23T14:49:24Z
- **Completed:** 2026-01-23T14:56:01Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Removed isicle_version field from JobStatus model with backwards compatibility
- Simplified create_job_directory() signature by removing isicle_version parameter
- Updated all API routes (jobs.py, web.py) to not pass isicle_version
- Old job status.json files with isicle_version can still load without errors

## Task Commits

Each task was committed atomically:

1. **Task 1: Remove isicle_version from JobStatus model** - `8b701c8` (refactor)
2. **Task 2: Update storage.py and API routers** - `b7e161c` (refactor)

## Files Created/Modified

- `src/qm_nmr_calc/models.py` - Added extra="ignore" to model_config, removed isicle_version field
- `src/qm_nmr_calc/storage.py` - Removed isicle_version parameter from create_job_directory()
- `src/qm_nmr_calc/api/routers/jobs.py` - Removed isicle_version="N/A" from job creation calls
- `src/qm_nmr_calc/api/routers/web.py` - Removed isicle_version="N/A" from job creation call

## Decisions Made

1. **Backwards compatibility strategy:** Used Pydantic's `extra="ignore"` in model_config to allow old job files with isicle_version field to load successfully while removing the field from the model definition
2. **Complete removal:** Removed isicle_version from all code paths (not just marked as deprecated) since ISiCLE is fully replaced by direct NWChem integration
3. **Comment update:** Changed "Versions (for reproducibility)" to "Version (for reproducibility)" (singular) to reflect only nwchem_version remains

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- ISiCLE version tracking completely removed from models and storage
- Only attribution comments remain in nwchem/__init__.py and runner.py
- API description in app.py still mentions "ISiCLE/NWChem" (could be updated to just "NWChem")
- Ready for final ISiCLE cleanup if needed (API description, remaining comments)
- All job storage and retrieval functions work correctly without isicle_version

**Note on test failures:** 5 tests fail due to pre-existing issue from 11-01 where `shielding_to_shift()` signature changed to require `functional`, `basis_set`, `solvent` parameters. These tests use the old signature and are unrelated to the isicle_version removal. Storage and model tests all pass.

---
*Phase: 11-production-integration*
*Completed: 2026-01-23*
