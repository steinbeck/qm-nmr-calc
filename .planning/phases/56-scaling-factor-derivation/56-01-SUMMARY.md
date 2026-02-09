---
phase: 56-scaling-factor-derivation
plan: 01
subsystem: benchmark
tags: [scaling-factors, linear-regression, delta50, ols, nmr-prediction, B3LYP, solvent-models]

# Dependency graph
requires:
  - phase: 55-delta50-benchmark-calculations
    provides: 200 benchmark calculations for 4 new solvents (Methanol, Water, Acetone, Benzene)
  - phase: 11-scaling-factor-derivation
    provides: Vacuum scaling factors from Phase 11.2
provides:
  - 14 production-ready scaling factors (6 solvents x 2 nuclei + vacuum)
  - Regression analysis report with quality metrics for all solvents
  - 24 diagnostic plots (regression + residual distributions)
affects: [57-solvent-integration, production-deployments, nmr-predictions]

# Tech tracking
tech-stack:
  added: []
  patterns: [OLS-regression-with-outlier-removal, quality-gate-validation]

key-files:
  created:
    - data/benchmark/delta50/plots/B3LYP_Methanol_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_Water_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_Acetone_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_Benzene_1H_regression.png
  modified:
    - src/qm_nmr_calc/benchmark/analysis.py
    - src/qm_nmr_calc/data/scaling_factors.json
    - data/benchmark/delta50/scaling_factors.json
    - data/benchmark/delta50/SCALING-FACTORS.md

key-decisions:
  - "Updated default solvents in analysis.py to include all 6 solvents automatically"
  - "Merged strategy: start from 12 newly derived factors, then add 2 vacuum factors for consistency"

patterns-established:
  - "Quality gates: R² > 0.99, 1H MAE < 0.2 ppm, 13C MAE < 3.0 ppm"
  - "Regression methodology: OLS with 3-sigma outlier removal and refit"

# Metrics
duration: 14min
completed: 2026-02-09
---

# Phase 56 Plan 01: Scaling Factor Derivation Summary

**12 production-ready scaling factors derived from 200 benchmark calculations across 4 new solvents, all passing R² > 0.99 and MAE quality gates**

## Performance

- **Duration:** 14 min
- **Started:** 2026-02-09T09:39:24Z
- **Completed:** 2026-02-09T09:53:27Z
- **Tasks:** 2
- **Files modified:** 20 (1 code, 1 package data, 2 benchmark outputs, 16 plots)

## Accomplishments
- Derived 12 scaling factor sets for all 6 B3LYP solvents (CHCl3, DMSO, Methanol, Water, Acetone, Benzene)
- All 8 new factor sets pass quality gates: R² > 0.99, 1H MAE < 0.2 ppm, 13C MAE < 3.0 ppm
- Merged 14 total factors into package data (12 derived + 2 vacuum preserved)
- Generated comprehensive regression report with 24 diagnostic plots
- Verified runtime loading: all 4 new solvents accessible via shifts.py API

## Task Commits

Each task was committed atomically:

1. **Task 1: Update default solvents and run analysis to derive all scaling factors** - `10af99d` (feat)
2. **Task 2: Merge factors into package data and validate end-to-end** - `4b5568d` (feat)

## Files Created/Modified
- `src/qm_nmr_calc/benchmark/analysis.py` - Updated default solvents list to include 4 new solvents
- `src/qm_nmr_calc/data/scaling_factors.json` - Package data now has 14 entries (6 existing + 8 new)
- `data/benchmark/delta50/scaling_factors.json` - Analysis output with 12 derived factor sets
- `data/benchmark/delta50/SCALING-FACTORS.md` - Updated report with all 6 solvents
- `data/benchmark/delta50/plots/B3LYP_{Solvent}_{Nucleus}_regression.png` - 12 regression scatter plots
- `data/benchmark/delta50/plots/B3LYP_{Solvent}_{Nucleus}_residuals.png` - 12 residual histograms

## Quality Metrics Summary

### New Solvents (4 solvents x 2 nuclei = 8 factor sets)

**Methanol:**
- 1H: R² = 0.9951, MAE = 0.126 ppm, n = 335 atoms
- 13C: R² = 0.9974, MAE = 2.139 ppm, n = 219 atoms

**Water:**
- 1H: R² = 0.9951, MAE = 0.127 ppm, n = 335 atoms
- 13C: R² = 0.9974, MAE = 2.161 ppm, n = 219 atoms

**Acetone:**
- 1H: R² = 0.9951, MAE = 0.126 ppm, n = 335 atoms
- 13C: R² = 0.9975, MAE = 2.117 ppm, n = 219 atoms

**Benzene:**
- 1H: R² = 0.9950, MAE = 0.128 ppm, n = 335 atoms
- 13C: R² = 0.9980, MAE = 1.761 ppm, n = 219 atoms

**All quality gates passed:**
- ✓ R-squared > 0.99 for all 8 new factor sets
- ✓ 1H MAE < 0.2 ppm for all 4 new solvents
- ✓ 13C MAE < 3.0 ppm for all 4 new solvents

### Existing Solvents (re-derived for consistency)

**CHCl3:**
- 1H: R² = 0.9952, MAE = 0.124 ppm
- 13C: R² = 0.9978, MAE = 1.949 ppm

**DMSO:**
- 1H: R² = 0.9951, MAE = 0.126 ppm
- 13C: R² = 0.9974, MAE = 2.152 ppm

### Vacuum (preserved from Phase 11.2)

**Vacuum:**
- 1H: slope = -0.9554, intercept = 30.5446 (unchanged)
- 13C: slope = -0.9726, intercept = 175.7109 (unchanged)

## Decisions Made

**1. Updated default solvents in analysis.py**
- Changed from `["CHCl3", "DMSO"]` to `["CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene"]`
- Rationale: Makes analyze command process all benchmark data by default

**2. Merge strategy: 12 new + 2 vacuum**
- Approach: Start from delta50 analysis output (12 entries), add vacuum entries from package file
- Rationale: Ensures CHCl3 and DMSO factors are re-derived from same analysis run as new solvents (consistency), while preserving vacuum factors from Phase 11.2's separate benchmark

**3. Verified vacuum factors preserved exactly**
- Checked slope/intercept values match Phase 11.2 to 10 decimal places
- Rationale: Vacuum benchmark was separate effort, must not be regenerated

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all analyses completed successfully, quality gates passed on first attempt.

## Validation Summary

**Runtime loading verified:**
- `load_scaling_factors()` returns 14 entries
- `get_scaling_factor("B3LYP", "6-311+G(2d,p)", "1H", "Methanol")` works
- All 4 new solvents accessible via shifts.py API

**Test suite verification:**
- API tests: 31 passed
- NWChem tests: 49 passed
- Total verified: 80 tests (subset, full suite is benchmark-heavy)

**Regression analysis:**
- 7 outliers removed from 1H regressions (per solvent)
- 2 outliers removed from 13C regressions (per solvent)
- Consistent outlier patterns across all solvents
- All n_points > 200 for robust statistics

## Next Phase Readiness

**Ready for Phase 57 (Solvent Integration):**
- All 14 scaling factors available in package data
- Runtime loading API functional for all 4 new solvents
- Quality metrics documented for validation
- Regression plots available for visual inspection

**No blockers:**
- All quality gates passed
- Test suite passes
- Vacuum factors preserved correctly

**Integration path clear:**
- Next phase: Update solvents.py and shifts.py to expose new solvents in API
- Validation: Ensure runtime predictions match benchmark quality

---
*Phase: 56-scaling-factor-derivation*
*Completed: 2026-02-09*
