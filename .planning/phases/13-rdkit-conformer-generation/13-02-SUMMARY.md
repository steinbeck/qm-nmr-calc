---
phase: 13-rdkit-conformer-generation
plan: "02"
subsystem: conformer-generation
tags: [rdkit, rmsd, conformer-filtering, energy-window, symmetry-aware]

# Dependency graph
requires:
  - phase: 12-conformer-data-model
    provides: ConformerData and ConformerEnsemble models for multi-conformer support
provides:
  - RMSD-based conformer deduplication using symmetry-aware comparison
  - Energy window filtering for conformer selection
  - Pure, reusable filter functions decoupled from RDKit Mol objects (energy filtering)
affects: [13-03-pipeline, 14-boltzmann, conformer-ensemble-workflow]

# Tech tracking
tech-stack:
  added: [rdkit.Chem.rdMolAlign.GetBestRMS]
  patterns: ["Greedy RMSD deduplication algorithm", "Decoupled energy filtering (operates on extracted data, not Mol objects)"]

key-files:
  created:
    - src/qm_nmr_calc/conformers/filters.py
    - tests/test_conformer_filters.py
  modified: []

key-decisions:
  - "Use GetBestRMS (pairwise) rather than GetAllConformerBestRMS (full matrix) - clearer code, marginally slower but acceptable for typical ensemble sizes"
  - "Decouple energy filtering from RDKit Mol objects - makes function reusable for MMFF and DFT energies with synthetic test data"
  - "Greedy algorithm for RMSD deduplication - keep first, skip similar to any kept conformer"

patterns-established:
  - "Conformer filtering functions are pure (no side effects) and deterministic"
  - "Energy filtering operates on parallel conf_ids and energies lists"
  - "RMSD comparison uses symmetry-aware GetBestRMS for accurate duplicate detection"

# Metrics
duration: 4min
completed: 2026-01-27
---

# Phase 13 Plan 02: Conformer Filters Summary

**Symmetry-aware RMSD deduplication and energy window filtering using greedy algorithm and decoupled energy filtering**

## Performance

- **Duration:** 4 min
- **Started:** 2026-01-27T07:44:53Z
- **Completed:** 2026-01-27T07:48:40Z
- **Tasks:** 1 (TDD - produced 2 commits: test + feat)
- **Files modified:** 2 (1 created implementation, 1 created test)

## Accomplishments
- Implemented `deduplicate_by_rmsd` using GetBestRMS for symmetry-aware comparison
- Implemented `filter_by_energy_window` decoupled from RDKit objects for reusability
- Comprehensive test suite with 12 test cases covering edge cases
- All tests pass with no regressions in existing test suite

## Task Commits

Each task was committed atomically following TDD RED-GREEN-REFACTOR:

1. **Task 1 (RED): Write failing tests** - `a71ec94` (test)
2. **Task 1 (GREEN): Implement filter functions** - `fea33c3` (feat)

**Plan metadata:** (to be committed after SUMMARY)

_Note: TDD task produced 2 commits (test → feat). No refactor phase needed - code was clean on first pass._

## Files Created/Modified
- `src/qm_nmr_calc/conformers/filters.py` - RMSD deduplication and energy window filtering functions
- `tests/test_conformer_filters.py` - Comprehensive test suite with 12 test cases

## Decisions Made

**Use pairwise GetBestRMS instead of matrix approach**
- GetAllConformerBestRMS returns flattened lower-triangle matrix requiring careful index math
- Pairwise GetBestRMS is clearer and only marginally slower for typical ensemble sizes (< 200 conformers)
- Prioritized code clarity for maintainability

**Decouple energy filtering from RDKit Mol objects**
- filter_by_energy_window accepts pre-extracted conf_ids and energies lists
- Makes function reusable for both MMFF and DFT energy filtering
- Enables testing with synthetic data without creating RDKit Mol objects
- Follows single-responsibility principle

**Greedy deduplication algorithm**
- Keep first conformer, compare each subsequent to all kept conformers
- If RMSD below threshold to ANY kept conformer, mark as duplicate
- Simple, deterministic, and effective for conformer deduplication

## Deviations from Plan

None - plan executed exactly as written.

TDD process followed: RED (failing tests) → GREEN (implementation) → verification (all tests pass).

## Issues Encountered

None - implementation was straightforward. GetBestRMS handles symmetry automatically, and energy filtering logic is simple arithmetic.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Phase 13 Plan 03 (Pipeline integration):**
- Filter functions tested and verified
- Both functions have clear APIs with docstrings and examples
- Pure functions with no side effects - easy to integrate into pipeline
- Energy filtering reusable for both MMFF (Plan 03) and DFT (future phases)

**Notes for next plan:**
- deduplicate_by_rmsd returns list of conformer IDs to keep
- Caller is responsible for removing filtered conformers from Mol object
- filter_by_energy_window preserves input order for deterministic results

---
*Phase: 13-rdkit-conformer-generation*
*Completed: 2026-01-27*
