---
phase: 64-solvent-integration
plan: 01
subsystem: core
tags: [solvents, nmr, scaling-factors, nmredata]

# Dependency graph
requires:
  - phase: 63-scaling-factor-derivation
    provides: Scaling factors for 6 new solvents (pyridine, thf, toluene, dcm, acetonitrile, dmf)
provides:
  - 6 new solvents integrated into core modules (solvents.py, shifts.py, nmredata.py)
  - Validation and display name support for new solvents
  - NMReData export support for new solvents
  - Test coverage for new solvents (8 new tests)
affects: [nmr-prediction, api, cli, ui-solvent-selector]

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/solvents.py
    - src/qm_nmr_calc/shifts.py
    - src/qm_nmr_calc/nmredata.py
    - tests/test_nmredata.py

key-decisions:
  - "Maintain consistent solvent naming across all three modules"
  - "Add tests for each new solvent mapping"

patterns-established:
  - "New solvents require updates to 3 files: solvents.py, shifts.py, nmredata.py"

# Metrics
duration: 9min
completed: 2026-02-11
---

# Phase 64 Plan 01: Solvent Integration Summary

**6 new solvents (pyridine, thf, toluene, dcm, acetonitrile, dmf) integrated into core modules with validation, scaling factor mapping, and NMReData export support**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-11T08:54:18Z
- **Completed:** 2026-02-11T09:03:35Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Expanded SUPPORTED_SOLVENTS from 7 to 13 entries
- Added scaling factor mappings for all 6 new solvents
- Added NMReData deuterated form mappings for all 6 new solvents
- Verified all 26 scaling factor lookups succeed (13 solvents x 2 nuclei)
- Added 6 new test methods for solvent mappings
- Fixed existing test to use ethanol instead of toluene
- All 456 tests passing

## Task Commits

Each task was committed atomically:

1. **Task 1: Add 6 new solvents to solvents.py and shifts.py** - `5f345b3` (feat)
2. **Task 2: Add 6 new solvents to nmredata.py and update tests** - `276532d` (feat)

## Files Created/Modified
- `src/qm_nmr_calc/solvents.py` - Added 6 entries to SUPPORTED_SOLVENTS dict
- `src/qm_nmr_calc/shifts.py` - Added 6 entries to solvent_map in get_scaling_factor()
- `src/qm_nmr_calc/nmredata.py` - Added 6 entries to solvent_map in map_solvent_to_nmredata()
- `tests/test_nmredata.py` - Added 6 new test methods and fixed test_unknown_solvent_raises_error

## Decisions Made
None - followed plan as specified

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All 6 new solvents fully integrated and tested
- Ready for UI updates to expose new solvents in dropdowns/selectors
- Ready for documentation updates to list expanded solvent support
- Ready for CLI updates to include new solvents in help text

---
*Phase: 64-solvent-integration*
*Completed: 2026-02-11*
