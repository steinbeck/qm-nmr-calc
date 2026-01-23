---
phase: 10-scaling-factors
plan: 01
subsystem: benchmark
tags: [statsmodels, ols-regression, pydantic, pandas, nmr-scaling]

# Dependency graph
requires:
  - phase: 09-benchmark-calculations
    provides: DELTA50 benchmark results with shielding data in shifts.json
provides:
  - ScalingFactor Pydantic model for regression metadata
  - RegressionData model for shielding-shift pairs
  - aggregate_regression_data() for data collection
  - fit_scaling_factors() for OLS regression with outlier removal
  - derive_all_factors() for automated factor derivation
affects: [10-02-PLAN, 11-production-integration]

# Tech tracking
tech-stack:
  added: []  # statsmodels, pandas already in project
  patterns:
    - "OLS regression with residual-based outlier removal (3 sigma threshold)"
    - "Composite factor keys: functional/basis_set/nucleus/solvent"

key-files:
  created:
    - src/qm_nmr_calc/benchmark/analysis.py
  modified:
    - src/qm_nmr_calc/benchmark/models.py
    - src/qm_nmr_calc/benchmark/__init__.py

key-decisions:
  - "numpy array indexing for statsmodels conf_int() results"
  - "Default to B3LYP functional only (WP04 data incomplete)"

patterns-established:
  - "ScalingFactor model stores all regression metadata including CIs"
  - "get_factor_key() for consistent composite key generation"

# Metrics
duration: 3min
completed: 2026-01-23
---

# Phase 10 Plan 01: Scaling Factor Derivation Summary

**OLS regression analysis deriving B3LYP NMR scaling factors from DELTA50 benchmark data with 0.99+ R-squared**

## Performance

- **Duration:** 3 min
- **Started:** 2026-01-23T13:00:31Z
- **Completed:** 2026-01-23T13:03:38Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments

- ScalingFactor and RegressionData Pydantic models for structured regression data
- aggregate_regression_data() pairs 563 shielding-shift values across DELTA50
- fit_scaling_factors() with outlier removal achieving R^2 > 0.99 for all combinations
- derive_all_factors() producing 4 factor sets (B3LYP x 2 solvents x 2 nuclei)

## Derived Scaling Factors

| Key | Slope | Intercept | R^2 | MAE (ppm) | n_points |
|-----|-------|-----------|-----|-----------|----------|
| B3LYP/6-311+G(2d,p)/1H/CHCl3 | -0.9375 | 29.92 | 0.9952 | 0.124 | 335 |
| B3LYP/6-311+G(2d,p)/13C/CHCl3 | -0.9497 | 172.69 | 0.9978 | 1.949 | 219 |
| B3LYP/6-311+G(2d,p)/1H/DMSO | -0.9323 | 29.73 | 0.9951 | 0.126 | 335 |
| B3LYP/6-311+G(2d,p)/13C/DMSO | -0.9429 | 171.77 | 0.9974 | 2.152 | 219 |

## Task Commits

Each task was committed atomically:

1. **Task 1: Add ScalingFactor model** - `a999b1e` (feat)
2. **Task 2: Create analysis.py with regression** - `22a2e53` (feat)
3. **Task 3: Test and fix CI extraction** - `f092872` (fix)

## Files Created/Modified

- `src/qm_nmr_calc/benchmark/models.py` - Added ScalingFactor and RegressionData models
- `src/qm_nmr_calc/benchmark/analysis.py` - Created with 4 functions for factor derivation
- `src/qm_nmr_calc/benchmark/__init__.py` - Export new models and analysis functions

## Decisions Made

- Used numpy array indexing for statsmodels conf_int() which returns ndarray not DataFrame
- Default to B3LYP functional only since WP04 benchmark calculations may be incomplete

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed statsmodels conf_int() extraction**
- **Found during:** Task 3 (end-to-end testing)
- **Issue:** Code used DataFrame iloc but conf_int() returns numpy array
- **Fix:** Changed ci.iloc[x, y] to ci[x, y] array indexing
- **Files modified:** src/qm_nmr_calc/benchmark/analysis.py
- **Verification:** derive_all_factors() runs successfully, returns 4 factors
- **Committed in:** f092872

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Minor implementation bug fixed during testing. No scope creep.

## Issues Encountered

None - implementation straightforward following RESEARCH.md patterns.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Analysis module complete with all required functions
- Ready for Phase 10 Plan 2: Report generation with plots and tables
- Factor data structures ready for publication-quality output

---
*Phase: 10-scaling-factors*
*Completed: 2026-01-23*
