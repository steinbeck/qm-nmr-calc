---
phase: 09-benchmark-calculations
plan: 02
subsystem: benchmark
tags: [nwchem-execution, nmr-calculations, delta50, b3lyp, noautoz]

# Dependency graph
requires:
  - phase: 09-benchmark-calculations-01
    provides: benchmark runner with status tracking and headless mode
provides:
  - 100 B3LYP NMR shift calculations (50 molecules x 2 solvents)
  - shifts.json files for all DELTA50 molecules
  - Calculated 1H and 13C chemical shifts
  - BENCHMARK-RESULTS.md summary report
affects: [10-scaling-factors, accuracy-analysis]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - noautoz directive for linear molecule geometry handling
    - Parallel NWChem execution (8 processes)

key-files:
  created:
    - .planning/phases/09-benchmark-calculations/BENCHMARK-RESULTS.md
  modified:
    - src/qm_nmr_calc/nwchem/input_gen.py
    - data/benchmark/results/*/B3LYP_*/shifts.json (100 files)
    - data/benchmark/results/status.json

key-decisions:
  - "B3LYP functional only (WP04 deferred)"
  - "noautoz as optional parameter (not default) for linear molecule fallback"
  - "8 parallel processes for throughput"

patterns-established:
  - "Optional noautoz parameter for AUTOZ failure recovery (not default)"

# Metrics
duration: ~2.5h total
completed: 2026-01-22
---

# Phase 9 Plan 02: DELTA50 Benchmark Execution Summary

**Executed 100 B3LYP NMR calculations across 50 DELTA50 molecules in 2 solvents**

## Performance

- **Duration:** ~2.5 hours total (including compound_14 re-run)
- **Started:** 2026-01-22T14:47:34Z
- **Completed:** 2026-01-22T~17:15:00Z
- **Calculations:** 100 (50 molecules x 2 solvents)
- **Average calc time:** ~109 seconds

## Accomplishments

- Pilot run validated (5 molecules, 20 calculations)
- Full headless benchmark executed for remaining molecules
- compound_14 (2-butyne) issue diagnosed and fixed with noautoz
- 100/100 shifts.json files generated (100% success rate)
- BENCHMARK-RESULTS.md summary report created

## Benchmark Results

### Coverage
| Functional | Solvent | Completed |
|------------|---------|-----------|
| B3LYP | CHCl3 | 50/50 |
| B3LYP | DMSO | 50/50 |
| WP04 | CHCl3 | 0/50 |
| WP04 | DMSO | 0/50 |

### Shift Statistics
- **1H shifts:** 684 values, range -1.50 to 9.27 ppm, avg 1.91 ppm
- **13C shifts:** 442 values, range -10.72 to 239.11 ppm, avg 84.46 ppm

## Issues Encountered

### compound_14 AUTOZ Failure
- **Problem:** NWChem's AUTOZ failed for 2-butyne (linear alkyne)
- **Error:** "insufficient internal variables: expected 23 got 24"
- **Solution:** Added `noautoz` directive to geometry blocks in input_gen.py
- **Resolution:** compound_14 completed successfully after fix

## Commits

1. **noautoz fix for input_gen.py** - Modified geometry optimization and shielding templates

## Files Created/Modified

- `src/qm_nmr_calc/nwchem/input_gen.py` - Added noautoz directive to geometry blocks
- `data/benchmark/results/compound_*/B3LYP_*/shifts.json` - 100 result files
- `.planning/phases/09-benchmark-calculations/BENCHMARK-RESULTS.md` - Summary report

## Deviations from Plan

- **WP04 functional not executed:** Only B3LYP was run (100 calculations instead of planned 200)
- **Manual intervention required:** compound_14 required code fix and manual re-run

## User Setup Required

None - all calculations executed automatically.

## Next Phase Readiness

- 100 B3LYP shifts ready for Phase 10 scaling factor derivation
- Experimental shifts available in DELTA50 dataset
- Linear regression analysis can proceed
- WP04 can be run later if comparative analysis needed

## Verification Checklist

- [x] Pilot run completed (5 molecules, 20 calculations)
- [x] Full run completed for B3LYP functional
- [x] status.json shows final state
- [x] BENCHMARK-RESULTS.md created
- [x] Results ready for Phase 10

---
*Phase: 09-benchmark-calculations*
*Completed: 2026-01-22*
