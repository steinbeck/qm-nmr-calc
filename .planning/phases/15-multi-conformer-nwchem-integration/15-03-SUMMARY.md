---
phase: 15-multi-conformer-nwchem-integration
plan: 03
subsystem: nwchem
tags: [nwchem, nmr, dft, ensemble, boltzmann, shielding]

# Dependency graph
requires:
  - phase: 15-02
    provides: run_conformer_dft_optimization and apply_post_dft_filter functions
  - phase: 14-02
    provides: average_ensemble_nmr for Boltzmann-weighted averaging
  - phase: 11-01
    provides: shielding_to_shift with scaling factors
provides:
  - run_conformer_nmr_calculations: NMR shielding loop for optimized conformers
  - run_ensemble_dft_and_nmr: Full DFT -> filter -> NMR orchestrator
affects: [16-crest-integration, 17-api-ensemble-mode]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Sequential NMR processing with per-conformer scratch isolation"
    - "No minimum success threshold for NMR (unlike DFT's 50%)"
    - "nmr_basis_set substitution in preset for NMR step"

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/nwchem/runner.py
    - src/qm_nmr_calc/nwchem/__init__.py
    - tests/test_conformer_nwchem.py

key-decisions:
  - "No minimum success threshold for NMR step (any successful results usable)"
  - "nmr_basis_set replaces basis_set in preset for NMR calculations"
  - "Ensemble conformers mutated in place (no reconstruction needed)"

patterns-established:
  - "NMR step uses skip_optimization=True with optimized geometry"
  - "Status transitions: optimized -> nmr_running -> nmr_complete"
  - "Pipeline returns (ensemble, nmr_results) for Boltzmann averaging"

# Metrics
duration: 9min
completed: 2026-01-27
---

# Phase 15 Plan 03: NMR Loop and Full Ensemble Orchestrator Summary

**Per-conformer NMR shielding loop with skip_optimization=True and full DFT->filter->NMR->Boltzmann orchestrator using existing scaling factors**

## Performance

- **Duration:** 9 min
- **Started:** 2026-01-27T21:28:59Z
- **Completed:** 2026-01-27T21:38:00Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- NMR shielding calculations run on each post-DFT-filtered conformer with skip_optimization=True
- Shielding-to-shift conversion via existing shielding_to_shift with scaling factors
- Full ensemble orchestrator ties together: DFT opt -> post-DFT filter -> NMR -> returns results for Boltzmann averaging
- Partial NMR failures handled gracefully with no minimum success threshold (unlike DFT's 50%)

## Task Commits

Each task was committed atomically:

1. **Task 1: Write failing tests for NMR loop and orchestrator** - `2d96abe` (test)
2. **Task 2: Implement run_conformer_nmr_calculations and run_ensemble_dft_and_nmr** - `836939f` (feat)

_TDD execution: RED (9 tests fail with ImportError) → GREEN (all 19 tests pass)_

## Files Created/Modified
- `src/qm_nmr_calc/nwchem/runner.py` - Added run_conformer_nmr_calculations and run_ensemble_dft_and_nmr functions
- `src/qm_nmr_calc/nwchem/__init__.py` - Export new functions
- `tests/test_conformer_nwchem.py` - 9 new test cases (19 total tests now)

## Decisions Made
- **No minimum success threshold for NMR:** Unlike DFT optimization (requires >50% success), any successful NMR results are usable since DFT already caught systematic failures
- **nmr_basis_set substitution:** NMR calculations use `preset["nmr_basis_set"]` instead of `preset["basis_set"]` for higher accuracy
- **In-place ensemble mutation:** Ensemble conformers are mutated during loops - no reconstruction needed since ConformerData objects are updated by reference

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed mock path for shielding_to_shift in tests**
- **Found during:** Task 2 (running tests after implementation)
- **Issue:** Tests patched `qm_nmr_calc.nwchem.runner.shielding_to_shift` but function is imported from `shifts.py`
- **Fix:** Changed mock path to `qm_nmr_calc.shifts.shielding_to_shift` (3 occurrences)
- **Files modified:** tests/test_conformer_nwchem.py
- **Verification:** All 19 tests pass
- **Committed in:** 836939f (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Minor test fix for correct mocking. No changes to implementation logic.

## Issues Encountered
None - implementation followed plan exactly as specified.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
Phase 15 complete! Multi-conformer NWChem integration pipeline fully implemented:
- ✅ Phase 15-01: DFT energy extraction (hartree native units)
- ✅ Phase 15-02: DFT optimization loop with post-DFT filtering
- ✅ Phase 15-03: NMR shielding loop and full orchestrator

**Ready for Phase 16 (CREST integration):** CREST will replace RDKit KDG as high-accuracy conformer generator. The run_ensemble_dft_and_nmr orchestrator is ready to consume conformer ensembles from any source.

**Blockers:** None
**Concerns:** CREST timeout handling (macrocycles can hang) - Phase 16 will implement timeout with RDKit fallback

---
*Phase: 15-multi-conformer-nwchem-integration*
*Completed: 2026-01-27*
