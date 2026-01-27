---
phase: 14-boltzmann-averaging
plan: "01"
subsystem: conformers
tags: [boltzmann-weights, statistical-mechanics, numerical-stability, tdd]

# Dependency graph
requires:
  - phase: 12-conformer-data-model
    provides: EnergyUnit type and ConformerData model with energy tracking
provides:
  - calculate_boltzmann_weights function with numerical stability via exp-normalize trick
  - Energy unit conversion (kcal_mol, hartree, kj_mol)
  - Temperature-dependent weight calculation
affects: [15-conformer-nmr-integration, ensemble-averaging]

# Tech tracking
tech-stack:
  added: []
  patterns: [exp-normalize-trick, statistical-mechanics-calculations]

key-files:
  created:
    - src/qm_nmr_calc/conformers/boltzmann.py
    - tests/test_boltzmann.py
  modified: []

key-decisions:
  - "exp-normalize trick (subtract min energy before exp) for numerical stability"
  - "Pure Python math.exp (no numpy) to maintain minimal dependencies"
  - "RT = 0.001987204 * temperature_k for gas constant in kcal/(mol*K)"

patterns-established:
  - "Boltzmann weights always sum to 1.0 (within 1e-10 tolerance)"
  - "Single conformer returns [1.0] shortcut"
  - "Relative energy differences matter, not absolute values"

# Metrics
duration: 6min
completed: 2026-01-27
---

# Phase 14 Plan 01: Boltzmann Averaging Summary

**Numerically stable Boltzmann weight calculation from DFT energies using exp-normalize trick, supporting three energy units (kcal_mol, hartree, kj_mol) with comprehensive test coverage**

## Performance

- **Duration:** 6 minutes
- **Started:** 2026-01-27T15:10:08Z
- **Completed:** 2026-01-27T15:17:07Z
- **Tasks:** 1 (TDD: RED + GREEN phases)
- **Files modified:** 2

## Accomplishments
- Implemented numerically stable Boltzmann weight calculation using exp-normalize trick
- Supports three energy units with proper conversion factors
- Comprehensive test suite with 20 test cases covering edge cases, numerical stability, and correctness
- Zero regressions in existing conformer tests

## Task Commits

Each task was committed atomically:

1. **Task 1 (RED): Write failing tests** - `8bef3fa` (test)
2. **Task 1 (GREEN): Implement Boltzmann calculation** - `48d3b8b` (feat)

_Note: TDD tasks have multiple commits (test → feat). No refactor phase needed - implementation was clean._

## Files Created/Modified
- `src/qm_nmr_calc/conformers/boltzmann.py` - Boltzmann weight calculation with energy conversion and numerical stability
- `tests/test_boltzmann.py` - 20 test cases covering basic behavior, energy conversion, numerical stability, temperature handling, edge cases, and return contract

## Decisions Made

**exp-normalize trick for numerical stability**
- Subtract minimum energy before exponentiation: `exp(-(E_i - E_min)/RT)`
- All exponent arguments <= 0 → prevents overflow
- Minimum energy conformer gets exp(0) = 1.0
- High energy conformers underflow to 0.0 gracefully

**Pure Python math.exp (no numpy)**
- Maintains minimal dependencies (project already uses pure Python + Pydantic)
- Sufficient performance for typical conformer counts (1-200 conformers)

**Physical constants**
- R_KCAL = 0.001987204 kcal/(mol*K) - Gas constant
- HARTREE_TO_KCAL = 627.5095 - Hartree to kcal/mol conversion
- KJ_TO_KCAL = 1.0/4.184 - kJ/mol to kcal/mol conversion

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - TDD workflow proceeded smoothly. All tests passed on first GREEN implementation.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Phase 15 (Conformer NMR Integration):**
- Boltzmann weight calculation function complete and tested
- Energy unit handling matches ConformerData.energy_unit field from Phase 12
- Numerical stability verified with extreme energy ranges
- Temperature parameter ready for ensemble.temperature_k integration

**Test coverage:**
- 20 test cases covering all specified behavior
- Edge cases: empty list, invalid temperature
- Numerical stability: extreme ranges, large absolute values
- Physical correctness: known case (0.59 kcal/mol → 73%/27% at 298.15K)
- Return contract: sum to 1.0, non-negative, monotonic

**No blockers.** Ready to integrate into conformer ensemble NMR averaging.

---
*Phase: 14-boltzmann-averaging*
*Completed: 2026-01-27*
