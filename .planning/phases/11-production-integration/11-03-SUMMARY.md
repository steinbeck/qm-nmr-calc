---
phase: 11-production-integration
plan: 03
subsystem: api
tags: [fastapi, pydantic, api-schema, metadata, branding, delta50]

# Dependency graph
requires:
  - phase: 10-scaling-factors
    provides: get_scaling_factor() for looking up DELTA50 metadata
  - phase: 11-01
    provides: DELTA50 scaling factors loaded at production runtime
  - phase: 11-02
    provides: ISiCLE version field removed from models
provides:
  - NMRResultsResponse schema with scaling factor metadata fields
  - API responses include DELTA50 source and expected MAE
  - Complete NWChem-only branding (no ISiCLE references)
affects: [api-documentation, user-transparency, future-accuracy-discussions]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - API metadata includes provenance (source) and expected accuracy (MAE)
    - Expected accuracy formatted as "+/- X.XX ppm" for user readability
    - Factor lookup uses uppercase functional (stored lowercase in job status)

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/api/schemas.py
    - src/qm_nmr_calc/api/routers/jobs.py
    - src/qm_nmr_calc/api/app.py
    - src/qm_nmr_calc/api/templates/base.html

key-decisions:
  - "Expected MAE as string format '+/- X.XX ppm' (not raw float) for user-friendly display"
  - "Functional uppercase conversion needed for factor lookup (stored lowercase)"
  - "No error handling for unsupported factor combinations (ValueError propagates)"

patterns-established:
  - "API transparency pattern: include calculation metadata, source, and expected accuracy"
  - "Factor lookup pattern: uppercase functional, use job_status.nmr_results parameters"

# Metrics
duration: 5min
completed: 2026-01-23
---

# Phase 11 Plan 03: API Metadata Summary

**API responses now include DELTA50 scaling factor metadata (source, expected MAE) and complete NWChem-only branding**

## Performance

- **Duration:** 5 min
- **Started:** 2026-01-23T14:58:56Z
- **Completed:** 2026-01-23T15:03:26Z
- **Tasks:** 3
- **Files modified:** 4

## Accomplishments
- API responses transparently show scaling factor source (DELTA50) and expected accuracy
- Expected MAE displayed in user-friendly "+/- X.XX ppm" format for both 1H and 13C
- ISiCLE branding removed from API description and web UI footer
- All user-facing elements now credit NWChem only

## Task Commits

Each task was committed atomically:

1. **Task 1: Add factor metadata to NMRResultsResponse schema** - `798c1b8` (feat)
2. **Task 2: Update jobs router to include factor metadata in response** - `e9fdc44` (feat)
3. **Task 3: Update API description and UI branding** - `2930e30` (refactor)

**Bug fix:** Update tests for signature change - `bd7d32d` (fix)

## Files Created/Modified
- `src/qm_nmr_calc/api/schemas.py` - Added scaling_factor_source, h1_expected_mae, c13_expected_mae fields
- `src/qm_nmr_calc/api/routers/jobs.py` - Import get_scaling_factor, populate metadata in both status and results endpoints
- `src/qm_nmr_calc/api/app.py` - Updated description from "via ISiCLE/NWChem" to "via NWChem"
- `src/qm_nmr_calc/api/templates/base.html` - Updated footer from "Powered by ISiCLE and NWChem" to "Powered by NWChem"

## Decisions Made

**1. Expected MAE as string format**
- Decision: Use "+/- X.XX ppm" string format (not raw float) for h1_expected_mae and c13_expected_mae
- Rationale: More user-friendly and explicitly shows units, matches how chemists report accuracy

**2. No error handling for unsupported factor combinations**
- Decision: Let ValueError propagate from get_scaling_factor() if combination not supported
- Rationale: Unsupported combinations indicate data integrity bugs in validation pipeline. Job submission already validates solvent and preset constrains functionals, so this should never happen in production. If it does, surfacing as 500 error is appropriate.

**3. Functional case conversion**
- Decision: Convert functional to uppercase before calling get_scaling_factor()
- Rationale: Job status stores functional lowercase (e.g., "b3lyp") but factor keys use uppercase (e.g., "B3LYP/6-311+G(2d,p)/1H/CHCl3")

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Updated tests for shielding_to_shift signature change**
- **Found during:** Task 3 (test suite verification)
- **Issue:** Tests calling shielding_to_shift() without required functional/basis_set/solvent parameters (signature changed in Phase 10)
- **Fix:** Added production defaults (B3LYP, 6-311+G(2d,p), CHCl3) to all test calls
- **Files modified:** tests/test_nwchem_integration.py, tests/test_nwchem_output.py
- **Verification:** All 5 previously failing tests now pass
- **Committed in:** bd7d32d

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Bug fix necessary for test suite to pass. No scope creep.

## Issues Encountered

None

## Next Phase Readiness

- API now provides full transparency about calculation accuracy
- Users can see expected error ranges for their results
- ISiCLE branding transition complete
- Ready for documentation phase (API docs can now showcase metadata features)

## Verification Summary

**Schema verification:**
- ✓ NMRResultsResponse includes scaling_factor_source field
- ✓ NMRResultsResponse includes h1_expected_mae field
- ✓ NMRResultsResponse includes c13_expected_mae field

**Response verification:**
- ✓ job_status_to_response() includes factor metadata
- ✓ get_nmr_results() endpoint includes factor metadata
- ✓ get_scaling_factor imported and used correctly
- ✓ Expected MAE formatted as "+/- X.XX ppm"

**Branding verification:**
- ✓ API description says "via NWChem" (no ISiCLE)
- ✓ Web UI footer says "Powered by NWChem" (no ISiCLE)
- ✓ Only ISiCLE references are attribution comments (models.py, nwchem/runner.py)

---
*Phase: 11-production-integration*
*Completed: 2026-01-23*
