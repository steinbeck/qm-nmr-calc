---
phase: 55-delta50-benchmark-calculations
plan: 02
subsystem: benchmark
tags: [nwchem, cosmo, acetone, benzene, nmr-shielding, delta50, b3lyp]

# Dependency graph
requires:
  - phase: 54-benchmark-infrastructure
    provides: "Benchmark CLI with Acetone/Benzene solvent support, COSMO parameters"
  - phase: 55-01
    provides: "Shared bug fixes (cwd isolation, direct keyword) validated in pilot"
provides:
  - "50 Acetone B3LYP NMR shielding calculations (shifts.json)"
  - "50 Benzene B3LYP NMR shielding calculations (shifts.json)"
  - "BENCHMARK-RESULTS-AB.md execution report"
  - "Phase 55 fully complete (all 200 calculations across 4 new solvents)"
affects: [56-scaling-factor-derivation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "NWChem cwd isolation: set subprocess cwd to scratch dir for per-calculation scratch file isolation"
    - "DFT direct keyword: required for stable CPHF convergence with COSMO solvation"

key-files:
  created:
    - ".planning/phases/55-delta50-benchmark-calculations/BENCHMARK-RESULTS-AB.md"
  modified:
    - "src/qm_nmr_calc/nwchem/runner.py"
    - "src/qm_nmr_calc/nwchem/input_gen.py"

key-decisions:
  - "Fix NWChem scratch file pollution by setting cwd to scratch directory"
  - "Add DFT direct keyword for CPHF convergence with COSMO"

patterns-established:
  - "NWChem cwd isolation: always run NWChem with cwd set to the input file's parent directory to prevent scratch file interference between sequential calculations"

# Metrics
duration: 12h 25m (mostly NWChem compute time, ~20 min active work)
completed: 2026-02-08
---

# Phase 55 Plan 02: Acetone & Benzene Benchmark Calculations Summary

**100 NWChem DELTA50 benchmark calculations (50 Acetone + 50 Benzene) completed with 0 failures, completing Phase 55 with all 200 calculations across 4 new solvents**

## Performance

- **Duration:** ~12h 25m wall clock (mostly unattended NWChem compute; ~20 min active work across 3 sessions)
- **Started:** 2026-02-08T07:50:02Z
- **Completed:** 2026-02-08T20:15:05Z
- **Tasks:** 3 (pilot run, headless launch, verification + report)
- **Files modified:** 3

## Accomplishments

- 50 Acetone B3LYP NMR shielding calculations completed (342 1H + 221 13C data points)
- 50 Benzene B3LYP NMR shielding calculations completed (342 1H + 221 13C data points)
- Fixed NWChem scratch file isolation bug that caused 80% pilot failure rate
- Fixed CPHF convergence issue with COSMO solvation via DFT direct keyword
- All 200 calculations across all 4 new solvents (Methanol, Water, Acetone, Benzene) verified complete -- Phase 55 fully satisfied

## Task Commits

Each task was committed atomically:

1. **Task 1: Pilot run for Acetone and Benzene** - `df3ca05` (fix) -- bug fixes + pilot validation
2. **Task 2: Launch full headless run** - (no code changes, execution-only task)
3. **Task 3: Verify completion and generate report** - `b918d78` (docs) -- BENCHMARK-RESULTS-AB.md

## Files Created/Modified

- `src/qm_nmr_calc/nwchem/runner.py` - Added cwd parameter to subprocess.call for scratch file isolation
- `src/qm_nmr_calc/nwchem/input_gen.py` - Added 'direct' keyword to DFT block for CPHF stability
- `.planning/phases/55-delta50-benchmark-calculations/BENCHMARK-RESULTS-AB.md` - Execution report with stats and quality verification

## Decisions Made

| Decision | Rationale |
|----------|-----------|
| Fix scratch file isolation via cwd (not NWChem scratch directive) | Simplest fix: set subprocess cwd to input file directory so molecule.* files stay isolated per-calculation. Alternative was NWChem scratch_dir directive but that requires changing input format. |
| Add DFT direct keyword globally for shielding input | Required for CPHF convergence with COSMO. Slight performance cost (recompute integrals) but eliminates disk-based integral read failures. Applied to all shielding calculations, not just benchmark. |

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] NWChem scratch file pollution causing "could not read mo vectors"**
- **Found during:** Task 1 (pilot run)
- **Issue:** Sequential NWChem calculations all wrote molecule.movecs to project root. Stale files from previous calculations caused read failures -- 8/10 pilot calculations failed.
- **Fix:** Added `cwd=scratch_cwd` parameter to `subprocess.call()` in `run_nwchem()`, isolating scratch files per-calculation directory
- **Files modified:** `src/qm_nmr_calc/nwchem/runner.py`
- **Verification:** Re-ran pilot -- 10/10 successful, no stale files in project root
- **Committed in:** df3ca05

**2. [Rule 1 - Bug] CPHF convergence failure with COSMO solvation**
- **Found during:** Task 1 (pilot run, pre-existing fix in working tree)
- **Issue:** NMR property calculations failed with `cphf_solve2: SCF residual greater than 1d-2` when using disk-based integrals with COSMO
- **Fix:** Added `direct` keyword to DFT block in shielding input generator to force on-the-fly integral evaluation
- **Files modified:** `src/qm_nmr_calc/nwchem/input_gen.py`
- **Verification:** All 100 Acetone/Benzene calculations completed successfully
- **Committed in:** df3ca05

---

**Total deviations:** 2 auto-fixed (both Rule 1 - Bug)
**Impact on plan:** Both fixes were essential for the benchmark to run at all. Without them, 80%+ of calculations would fail. No scope creep -- fixes are minimal and targeted.

## Issues Encountered

- Initial pilot run had 80% failure rate (8/10) due to two NWChem bugs. Root cause analysis identified scratch file pollution and CPHF convergence issues. Both fixed before proceeding with full run.
- The "Authorization required, but no authorization protocol specified" messages in NWChem stderr are X11 display warnings (harmless), not actual authentication failures.

## User Setup Required

None - no external service configuration required.

## Benchmark Statistics

| Metric | Acetone | Benzene |
|--------|---------|---------|
| Calculations | 50/50 | 50/50 |
| Avg calc time | 295.3s (~4.9 min) | 295.8s (~4.9 min) |
| Total compute | 4.2 hours | 4.2 hours |
| 1H data points | 342 | 342 |
| 13C data points | 221 | 221 |
| 1H shielding range | 21.4 to 31.6 ppm | 21.5 to 31.6 ppm |
| 13C shielding range | -56.1 to 185.7 ppm | -50.3 to 185.1 ppm |

## Next Phase Readiness

- All 200 benchmark calculations (4 solvents x 50 molecules) complete with shielding data
- Phase 56 can proceed to derive scaling factors via linear regression on shielding vs experimental shift data
- No blockers or concerns

---
*Phase: 55-delta50-benchmark-calculations*
*Completed: 2026-02-08*
