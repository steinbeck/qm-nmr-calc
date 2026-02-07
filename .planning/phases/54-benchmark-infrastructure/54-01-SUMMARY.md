---
phase: 54-benchmark-infrastructure
plan: 01
subsystem: benchmark
tags: [nwchem, cosmo, solvents, benchmark, cli]

# Dependency graph
requires:
  - phase: 08-delta50-setup
    provides: benchmark runner infrastructure with CHCl3/DMSO solvents
provides:
  - "Benzene added to NWChem COSMO SUPPORTED_SOLVENTS"
  - "Benchmark SOLVENTS expanded to 6: CHCl3, DMSO, Methanol, Water, Acetone, Benzene"
  - "CLI --solvents accepts all 6 solvents"
affects: [55-benchmark-calculations, 56-scaling-factor-analysis, 57-solvent-integration]

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/nwchem/input_gen.py
    - src/qm_nmr_calc/benchmark/runner.py
    - src/qm_nmr_calc/benchmark/__main__.py
    - tests/test_nwchem_input.py
    - tests/test_benchmark.py

key-decisions:
  - "Benzene is the only new solvent needing addition to input_gen SUPPORTED_SOLVENTS (water, acetone, methanol were already there)"
  - "Runner converts title-case CLI names to lowercase via solvent.lower() for NWChem COSMO compatibility"

patterns-established: []

# Metrics
duration: 12min
completed: 2026-02-07
---

# Phase 54 Plan 01: Expanded Solvent Support Summary

**Added benzene to NWChem COSMO, expanded benchmark runner and CLI to 6 solvents (CHCl3, DMSO, Methanol, Water, Acetone, Benzene) for Phase 55 calculations**

## Performance

- **Duration:** 12 min
- **Started:** 2026-02-07T20:19:30Z
- **Completed:** 2026-02-07T20:31:31Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Added benzene to NWChem input_gen SUPPORTED_SOLVENTS, enabling COSMO solvation for benzene
- Expanded benchmark SOLVENTS list from 2 to 6, adding Methanol, Water, Acetone, Benzene
- Updated CLI --solvents choices to accept all 6 solvents without argparse errors
- Added 9 new test cases covering benzene COSMO generation, expanded solvent lists, and task matrix

## Task Commits

Each task was committed atomically:

1. **Task 1: Add 4 new solvents to input_gen, benchmark runner, and CLI** - `6f87f73` (feat)
2. **Task 2: Add tests for new solvent support** - `669364e` (test)

## Files Created/Modified
- `src/qm_nmr_calc/nwchem/input_gen.py` - Added "benzene" to SUPPORTED_SOLVENTS set
- `src/qm_nmr_calc/benchmark/runner.py` - Expanded SOLVENTS list to 6 entries, cleaned up comment
- `src/qm_nmr_calc/benchmark/__main__.py` - Expanded CLI --solvents choices, updated help text
- `tests/test_nwchem_input.py` - Added TestBenzeneSolvent class with 9 tests (3 benzene-specific + 6 parametrized)
- `tests/test_benchmark.py` - Added TestExpandedSolvents class with 6 tests

## Decisions Made
- Only benzene needed to be added to input_gen.py SUPPORTED_SOLVENTS; water, acetone, and methanol were already supported at that layer but not reachable from the benchmark CLI
- Kept title-case naming convention in SOLVENTS list (matching existing CHCl3/DMSO pattern) with runner converting to lowercase for NWChem

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All 6 solvents now accepted at every layer: CLI, runner, and NWChem input generation
- Phase 55 can dispatch `python -m qm_nmr_calc.benchmark run --solvents Methanol Water Acetone Benzene --functionals B3LYP` without errors
- 399 passing tests (5 pre-existing failures in NWChem integration and conformer ensemble tests, unrelated to this change)

---
*Phase: 54-benchmark-infrastructure*
*Completed: 2026-02-07*
