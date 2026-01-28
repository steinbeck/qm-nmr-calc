---
phase: 16-crest-integration
plan: 03
subsystem: conformers
tags: [crest, xtb, rdkit, conformer-generation, dispatch, health-check]

# Dependency graph
requires:
  - phase: 16-01
    provides: CREST binary detection and ALPB solvent mapping
  - phase: 16-02
    provides: CREST subprocess runner and ensemble builder
provides:
  - Unified conformer pipeline with CREST/RDKit dispatch
  - Health endpoint with CREST availability reporting
  - Package-level exports for CREST functions
affects: [17-api-integration, future-conformer-enhancements]

# Tech tracking
tech-stack:
  added: []
  patterns: [method-based-dispatch, fail-fast-validation, health-capability-reporting]

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/conformers/pipeline.py
    - src/qm_nmr_calc/conformers/__init__.py
    - src/qm_nmr_calc/api/routers/health.py
    - tests/test_conformer_pipeline.py
    - tests/test_api.py

key-decisions:
  - "conformer_method parameter for dispatch - 'rdkit_kdg' (default) or 'crest'"
  - "Fail-fast validation for CREST availability and solvent support"
  - "CREST detection as informational field in health endpoint"

patterns-established:
  - "Method-based dispatch: Single pipeline entry point routes to CREST or RDKit based on conformer_method parameter"
  - "Fail-fast validation: Check CREST availability and solvent support before processing to provide clear error messages"
  - "Health endpoint capability reporting: Non-blocking checks for optional features"

# Metrics
duration: 13min
completed: 2026-01-28
---

# Phase 16 Plan 03: CREST Integration Summary

**Unified conformer pipeline dispatches to CREST or RDKit based on method parameter with fail-fast validation and health endpoint capability reporting**

## Performance

- **Duration:** 13 min
- **Started:** 2026-01-28T11:07:43Z
- **Completed:** 2026-01-28T11:21:05Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Pipeline dispatches to CREST when method='crest' and CREST available
- Fail-fast validation prevents silent failures (checks binary availability and solvent support)
- Health endpoint exposes CREST availability for client decision-making
- Package exports enable direct import of CREST functions
- Full backward compatibility with existing RDKit-only path

## Task Commits

Each task was committed atomically:

1. **Task 1: Add CREST dispatch to pipeline and health endpoint** - `012ff36` (feat)
2. **Task 2: Add dispatch and health tests** - `9bcee38` (test)

## Files Created/Modified
- `src/qm_nmr_calc/conformers/pipeline.py` - Added conformer_method parameter with CREST/RDKit dispatch logic and fail-fast validation
- `src/qm_nmr_calc/conformers/__init__.py` - Exported detect_crest_available and generate_crest_ensemble
- `src/qm_nmr_calc/api/routers/health.py` - Added crest_available field to readiness endpoint
- `tests/test_conformer_pipeline.py` - Added TestCRESTDispatch class with 5 tests for dispatch and validation
- `tests/test_api.py` - Added test for crest_available field in health endpoint

## Decisions Made

**1. conformer_method parameter for dispatch**
- Added `conformer_method: str = "rdkit_kdg"` parameter to `generate_conformer_ensemble()`
- Values: "rdkit_kdg" (default) or "crest"
- Rationale: Explicit method selection at call site, backward compatible default, clear user intent

**2. Fail-fast validation for CREST**
- Check `detect_crest_available()` immediately when method='crest'
- Check `get_alpb_solvent()` returns non-None before proceeding
- Raise ValueError with clear error messages
- Rationale: Catch configuration errors early with actionable error messages rather than cryptic subprocess failures

**3. CREST as informational health field**
- CREST availability is informational only, doesn't block readiness (status 503)
- Exposed as both `checks["crest_available"]` and top-level `crest_available` field
- Rationale: CREST is optional enhancement, not required for basic operation. Clients can check availability and adjust conformer_method accordingly.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## Next Phase Readiness

**Ready for Phase 17 (API Integration):**
- Conformer pipeline has stable interface with method dispatch
- Health endpoint reports CREST availability for client decision-making
- All validation happens at pipeline entry point (no surprise failures mid-processing)
- Full test coverage for dispatch logic (5 new tests)

**Integration points for Phase 17:**
- API should accept `conformer_method` parameter in job submission
- API should query health endpoint to determine available methods
- API should pass solvent parameter to pipeline (already available for DFT, now used by CREST)
- API should handle ValueError exceptions from fail-fast validation and return 422 to client

---
*Phase: 16-crest-integration*
*Completed: 2026-01-28*
