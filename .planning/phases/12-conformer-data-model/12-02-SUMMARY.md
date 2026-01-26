---
phase: 12-conformer-data-model
plan: 02
subsystem: core
tags: [rdkit, atom-ordering, canonicalization, nwchem, tdd]

# Dependency graph
requires:
  - phase: 07-nwchem-integration
    provides: NWChem input/output handling and geometry utilities
provides:
  - Canonical atom ordering functions for deterministic atom indexing
  - NWChem to canonical index mapping
  - Atom count validation utilities
affects: [13-rdkit-conformer-generation, 14-boltzmann-weighting, 15-conformer-nmr-workflow]

# Tech tracking
tech-stack:
  added: []
  patterns: [canonical-atom-ordering, tdd-red-green-refactor]

key-files:
  created:
    - src/qm_nmr_calc/atom_ordering.py
    - tests/test_atom_ordering.py
  modified: []

key-decisions:
  - "Use RDKit's CanonicalRankAtoms for deterministic ordering"
  - "Keep mapping as metadata, don't physically reorder atoms (breaks RDKit operations)"
  - "NWChem uses 1-based indexing, canonical ranks are 0-based"

patterns-established:
  - "TDD with RED-GREEN-REFACTOR: separate commits for test/feat/refactor"
  - "Canonical ordering as ground truth for all conformer data"

# Metrics
duration: 3min
completed: 2026-01-26
---

# Phase 12 Plan 02: Canonical Atom Ordering Summary

**RDKit canonical atom ordering with deterministic 1-based NWChem mapping, ensuring atom index consistency across conformer lifecycle**

## Performance

- **Duration:** 3 min
- **Started:** 2026-01-26T15:10:18Z
- **Completed:** 2026-01-26T15:13:17Z
- **Tasks:** 1 (TDD: 3 commits)
- **Files modified:** 2

## Accomplishments
- Implemented canonical_atom_order() using RDKit's CanonicalRankAtoms for deterministic ordering
- Built map_nwchem_to_canonical() to translate 1-based NWChem indices to canonical ranks
- Created get_atom_count() validation helper for lifecycle consistency checks
- Achieved 100% test coverage with 13 passing tests
- Verified determinism over 100 runs for multiple molecules
- Confirmed symmetric molecules (benzene) handled correctly

## Task Commits

Each task was committed atomically following TDD:

1. **Task 1: Canonical Atom Ordering** (TDD cycle)
   - `e2b4bbb` (test) - RED: failing tests for atom ordering
   - `50caa1f` (feat) - GREEN: implementation with all tests passing
   - `efe4ce9` (refactor) - REFACTOR: extracted common SMILES parsing helper

**Plan metadata:** [will be added in final commit]

_TDD methodology: RED (failing test) → GREEN (minimal implementation) → REFACTOR (cleanup)_

## Files Created/Modified
- `src/qm_nmr_calc/atom_ordering.py` - Canonical ordering functions using RDKit's CanonicalRankAtoms
- `tests/test_atom_ordering.py` - Comprehensive test suite (13 tests covering determinism, mapping, symmetric molecules)

## Decisions Made

**1. Use RDKit's CanonicalRankAtoms**
- Rationale: Built-in, deterministic, well-tested canonical ordering algorithm
- Alternative considered: Custom ordering logic (rejected - reinventing the wheel)

**2. Keep mapping as metadata, don't physically reorder**
- Rationale: RenumberAtoms() would change mol object and break downstream RDKit operations
- Impact: Mapping is a translation layer, not a transformation

**3. Extract _smiles_to_mol_with_h() helper**
- Rationale: DRY principle - three functions all parse SMILES and add hydrogens
- Impact: 14 lines reduced to single helper function call

## Deviations from Plan

None - plan executed exactly as written. TDD red-green-refactor cycle produced three commits as expected.

## Issues Encountered

None - RDKit's CanonicalRankAtoms worked as expected for all test cases including symmetric molecules.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Phase 12 Plan 03 (RDKit Conformer Generation):**
- Canonical ordering functions available for conformer atom index tracking
- get_atom_count() provides validation for atom count consistency
- map_nwchem_to_canonical() ready for NMR output parsing

**No blockers.** Atom ordering foundation is solid and tested.

**Key insight for next phases:**
- All conformer data should reference canonical atom indices as ground truth
- RDKit's atom order after AddHs is deterministic for a given SMILES
- NWChem output uses 1-based indexing - use map_nwchem_to_canonical() for translation

---
*Phase: 12-conformer-data-model*
*Completed: 2026-01-26*
