---
phase: 55-delta50-benchmark-calculations
plan: 01
subsystem: benchmark
tags: [nwchem, cosmo, methanol, water, b3lyp, nmr-shielding, delta50]

# Dependency graph
requires:
  - phase: 54-benchmark-infrastructure
    provides: Benchmark CLI with expanded solvent support (Methanol, Water, Acetone, Benzene)
provides:
  - 50 Methanol B3LYP NMR shielding calculations (342 H + 221 C data points)
  - 50 Water B3LYP NMR shielding calculations (342 H + 221 C data points)
  - BENCHMARK-RESULTS-MW.md execution summary
  - CPHF convergence fix (direct integral evaluation for NMR+COSMO)
affects: [56-scaling-factor-derivation, 55-02-benchmark-calculations-ab]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "NWChem shielding input requires 'direct' keyword for CPHF convergence with COSMO solvation"
    - "NWChem scratch file isolation via cwd per-calculation prevents movecs pollution"

key-files:
  created:
    - .planning/phases/55-delta50-benchmark-calculations/BENCHMARK-RESULTS-MW.md
    - data/benchmark/results/compound_*/B3LYP_Methanol/shifts.json (50 files)
    - data/benchmark/results/compound_*/B3LYP_Water/shifts.json (50 files)
  modified:
    - src/qm_nmr_calc/nwchem/input_gen.py
    - src/qm_nmr_calc/nwchem/runner.py

key-decisions:
  - "Added 'direct' keyword to DFT shielding input for CPHF convergence with COSMO"
  - "Set NWChem cwd to scratch directory for per-calculation file isolation"

patterns-established:
  - "NMR shielding + COSMO requires 'direct' in DFT block to avoid cphf_solve2 SCF residual errors"

# Metrics
duration: ~8.5 hours compute (managed across 3 sessions)
completed: 2026-02-08
---

# Phase 55 Plan 01: DELTA50 Methanol & Water Benchmarks Summary

**100 NWChem B3LYP/6-311+G(2d,p) NMR shielding calculations completed for Methanol and Water solvents across 50 DELTA50 molecules with 0 failures, after fixing CPHF convergence bug requiring 'direct' integral evaluation**

## Performance

- **Duration:** ~8.5 hours NWChem compute (managed across 3 checkpoint sessions)
- **Started:** 2026-02-08T08:50:00Z (pilot run)
- **Completed:** 2026-02-08T17:44:39Z (full run finished)
- **Tasks:** 3 (pilot run + headless launch + verification/summary)
- **Files modified:** 2 source files + 100 result files + 1 report

## Accomplishments

- Fixed critical CPHF convergence bug: added `direct` keyword to NMR shielding DFT block, enabling NMR property calculations with COSMO solvation
- Fixed scratch file isolation: NWChem cwd set to input directory preventing cross-calculation movecs pollution
- Completed 50/50 Methanol and 50/50 Water benchmark calculations with zero failures
- Validated shielding data quality: 342 H and 221 C data points per solvent in physically reasonable ranges
- BENCH-02 (Methanol) and BENCH-03 (Water) success criteria satisfied

## Task Commits

Each task was committed atomically:

1. **Task 1: Pilot run (5 molecules x Methanol + Water)** - `df3ca05` (fix: CPHF convergence + scratch isolation)
2. **Task 2: Launch full headless run** - (no code changes, runtime-only task)
3. **Task 3: Verify completion + generate BENCHMARK-RESULTS-MW.md** - `53b2414` (feat)

## Files Created/Modified

- `src/qm_nmr_calc/nwchem/input_gen.py` - Added `direct` keyword to DFT shielding block for CPHF convergence
- `src/qm_nmr_calc/nwchem/runner.py` - Set cwd to scratch directory for per-calculation file isolation
- `data/benchmark/results/compound_*/B3LYP_Methanol/shifts.json` - 50 Methanol shielding result files
- `data/benchmark/results/compound_*/B3LYP_Water/shifts.json` - 50 Water shielding result files
- `.planning/phases/55-delta50-benchmark-calculations/BENCHMARK-RESULTS-MW.md` - Execution summary report

## Decisions Made

| Decision | Rationale |
|----------|-----------|
| Added `direct` to DFT shielding input | Required for CPHF convergence with COSMO; without it all calcs fail with SCF residual > 1d-2 |
| Set NWChem cwd to input file parent dir | Prevents molecule.movecs and molecule.db cross-contamination between sequential calculations |
| Kept ISiCLE-compatible format (named solvents, do_gasphase False) | Maintains compatibility with optimization inputs; only added `direct` for shielding |

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] CPHF convergence failure with COSMO solvation**
- **Found during:** Task 1 (pilot run)
- **Issue:** All 10 pilot calculations failed with `cphf_solve2: SCF residual greater than 1d-2` (0.41 for Methanol, 0.33 for Water). NWChem's CPHF solver requires direct integral evaluation when used with COSMO solvation for NMR properties.
- **Fix:** Added `direct` keyword to DFT block in `generate_shielding_input()`. This forces on-the-fly integral computation rather than cached integrals, which is necessary for stable CPHF convergence with COSMO.
- **Files modified:** `src/qm_nmr_calc/nwchem/input_gen.py`
- **Verification:** Re-ran pilot: 10/10 calculations succeeded. Full run: 100/100 succeeded.
- **Committed in:** df3ca05

**2. [Rule 3 - Blocking] NWChem scratch file cross-contamination**
- **Found during:** Task 1 (pilot run debugging)
- **Issue:** Sequential NWChem runs in the same working directory share molecule.movecs and molecule.db files, causing "could not read mo vectors" errors on subsequent calculations.
- **Fix:** Set `cwd` parameter in `subprocess.call()` to the input file's parent directory, isolating scratch files per-calculation.
- **Files modified:** `src/qm_nmr_calc/nwchem/runner.py`
- **Verification:** All 100 sequential calculations completed without movecs errors.
- **Committed in:** df3ca05

---

**Total deviations:** 2 auto-fixed (1 bug, 1 blocking)
**Impact on plan:** Both fixes were essential for correct NWChem operation. Without the `direct` fix, zero calculations would succeed. Without scratch isolation, sequential runs would fail intermittently.

## Issues Encountered

- Initial pilot run had 100% failure rate due to missing `direct` keyword in NMR shielding DFT input. Root cause analysis revealed the keyword was present in the old input format (pre-ISiCLE migration) but was dropped during the format change. Fixed by adding it back.
- The `show_status()` display shows combined progress from multiple concurrent benchmark processes sharing the same status.json file, which can be confusing but does not affect results.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Methanol and Water shielding data ready for scaling factor derivation in Phase 56
- Plan 55-02 (Acetone + Benzene benchmarks) was running concurrently and also completed
- All 200 calculations across 4 solvents available for Phase 56 analysis
- No blockers

---
*Phase: 55-delta50-benchmark-calculations*
*Completed: 2026-02-08*
