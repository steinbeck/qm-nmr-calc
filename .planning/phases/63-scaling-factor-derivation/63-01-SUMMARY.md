---
phase: 63
plan: 01
subsystem: benchmark
tags: [scaling-factors, regression, analysis, ols, quality-gates]
dependencies:
  requires:
    - 60-01-PLAN.md  # DELTA50 Pyridine & THF benchmark data
    - 61-01-PLAN.md  # DELTA50 Toluene & DCM benchmark data
    - 62-01-PLAN.md  # DELTA50 Acetonitrile & DMF benchmark data
  provides:
    - 26 scaling factor sets (13 solvents x 2 nuclei)
    - Quality-validated factors (R² > 0.99, MAE within gates)
    - Updated SCALING-FACTORS.md report with plots
  affects:
    - 64-01-PLAN.md  # Solvent integration (will use these factors at runtime)
tech-stack:
  added: []
  patterns: []
key-files:
  created:
    - data/benchmark/delta50/plots/B3LYP_Pyridine_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_Pyridine_13C_regression.png
    - data/benchmark/delta50/plots/B3LYP_THF_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_THF_13C_regression.png
    - data/benchmark/delta50/plots/B3LYP_Toluene_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_Toluene_13C_regression.png
    - data/benchmark/delta50/plots/B3LYP_DCM_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_DCM_13C_regression.png
    - data/benchmark/delta50/plots/B3LYP_Acetonitrile_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_Acetonitrile_13C_regression.png
    - data/benchmark/delta50/plots/B3LYP_DMF_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_DMF_13C_regression.png
    - data/benchmark/delta50/plots/B3LYP_*_residuals.png
  modified:
    - src/qm_nmr_calc/benchmark/analysis.py
    - src/qm_nmr_calc/data/scaling_factors.json
    - data/benchmark/delta50/scaling_factors.json
    - data/benchmark/delta50/SCALING-FACTORS.md
decisions:
  MERGE-STRATEGY: "Start from newly generated analysis output (24 solvent entries), add vacuum entries from existing package data (2 entries), to ensure all solvent factors come from same analysis run (consistency)"
  QUALITY-GATES: "All 6 new solvents pass R² > 0.99, 1H MAE < 0.2 ppm, 13C MAE < 3.0 ppm"
metrics:
  duration: 18m 38s
  completed: 2026-02-11
  tasks: 2
  commits: 2
  tests: 448 passed, 2 skipped
---

# Phase 63 Plan 01: Scaling Factor Derivation Summary

**One-liner:** Derived OLS scaling factors for 6 new solvents (Pyridine, THF, Toluene, DCM, Acetonitrile, DMF), validated all quality gates, and merged into package data with 26 total factor sets.

## Context

Phases 60-62 produced 300 benchmark calculations (shielding values) across 6 new solvents. This phase converts those raw shielding values into usable scaling factors via linear regression. Scaling factors are the critical link between calculated NMR shielding (sigma) and experimental chemical shifts (delta) - without them, the benchmark data is useless for predictions.

The task: Run OLS regression for each solvent x nucleus combination (12 new factor sets), validate they meet quality standards, and integrate all 26 factor sets (14 existing + 12 new) into the package data that gets loaded at runtime.

## Tasks Completed

### Task 1: Update default solvents and run analysis to derive all scaling factors

**Goal:** Update analysis.py to include all 12 solvents, run the analysis command to derive factors and generate report.

**Implementation:**
- Updated `derive_all_factors()` default solvents parameter from 6 to 12 solvents:
  - Existing: CHCl3, DMSO, Methanol, Water, Acetone, Benzene
  - New: Pyridine, THF, Toluene, DCM, Acetonitrile, DMF
- Ran `python -m qm_nmr_calc.benchmark analyze --output-dir data/benchmark/delta50`
- Analysis performed OLS regression with outlier removal (3-sigma threshold)
- Generated 24 plots (regression + residuals for each factor set)

**Key findings:**
- All 24 factor sets derived successfully (12 solvents x 2 nuclei)
- Outlier removal: 7 outliers removed from each 1H regression, 2 from each 13C regression
- Consistent behavior across all 12 solvents (similar R², MAE, RMSD values)

**Outputs:**
- `data/benchmark/delta50/scaling_factors.json` (24 entries)
- `data/benchmark/delta50/SCALING-FACTORS.md` (comprehensive report)
- `data/benchmark/delta50/plots/` (48 plots: 24 regression + 24 residual)

**Commit:** `e5cd614` - feat(63-01): derive scaling factors for 6 new solvents

### Task 2: Merge factors into package data and validate end-to-end

**Goal:** Merge analysis output into package scaling_factors.json, preserving vacuum entries, and validate with test suite.

**Implementation:**
- Created merge script using orjson library
- Merge strategy: Start with 24 newly derived solvent factors, add 2 vacuum factors from existing package data
- Verified vacuum factors unchanged (exact match to 10+ decimal places)
- Validated all 26 entries have required fields
- Ran full test suite: 448 tests passed, 2 skipped, 0 failures

**Key validations:**
- Total entries: 26 (24 solvent + 2 vacuum)
- `load_scaling_factors()` returns all 26 entries
- All 6 new solvent keys accessible for both 1H and 13C
- No test regressions introduced

**Commit:** `f3e2384` - feat(63-01): merge all 26 scaling factors into package data

## Quality Gates Validation

All 6 new solvents passed quality gates:

| Solvent       | 1H R²   | 1H MAE (ppm) | 13C R²  | 13C MAE (ppm) |
|---------------|---------|--------------|---------|---------------|
| Pyridine      | 0.9952  | 0.125 ✓      | 0.9976  | 2.085 ✓       |
| THF           | 0.9952  | 0.124 ✓      | 0.9977  | 2.023 ✓       |
| Toluene       | 0.9951  | 0.128 ✓      | 0.9980  | 1.774 ✓       |
| DCM           | 0.9952  | 0.124 ✓      | 0.9976  | 2.048 ✓       |
| Acetonitrile  | 0.9951  | 0.126 ✓      | 0.9974  | 2.143 ✓       |
| DMF           | 0.9951  | 0.126 ✓      | 0.9974  | 2.144 ✓       |

**Quality gates:**
- R² > 0.99: ✓ All 12 factor sets (range: 0.9951-0.9980)
- 1H MAE < 0.2 ppm: ✓ All 6 solvents (range: 0.124-0.128 ppm)
- 13C MAE < 3.0 ppm: ✓ All 6 solvents (range: 1.774-2.144 ppm)

**Observations:**
- Aromatic solvents (Benzene, Toluene) show best 13C performance (MAE ~1.76 ppm)
- Polar solvents show consistent 1H performance (MAE 0.124-0.127 ppm)
- R² values tightly clustered (0.9951-0.9980), indicating excellent fit quality
- All 6 new solvents on par with existing 6 solvents (no degradation)

## Deviations from Plan

None - plan executed exactly as written.

## Decisions Made

**MERGE-STRATEGY:**
- **Decision:** Start from newly generated analysis output (24 entries), then add vacuum entries from existing package data (2 entries)
- **Rationale:** Ensures all solvent factors come from same analysis run (consistency), while preserving vacuum factors from Phase 11.2's separate benchmark
- **Alternative considered:** Merge incrementally (add only 6 new solvents to existing 14) - rejected because would mix factor sets from different analysis runs
- **Impact:** All 12 solvent factors re-derived together, ensuring statistical consistency

**QUALITY-GATES:**
- **Decision:** Strict validation (R² > 0.99, 1H MAE < 0.2 ppm, 13C MAE < 3.0 ppm)
- **Rationale:** These gates established in Phase 11.2 ensure factors suitable for production NMR prediction
- **Result:** All 6 new solvents pass all gates on first attempt
- **Impact:** No need for outlier investigation or benchmark recalculation

## Technical Notes

### OLS Regression Methodology

The analysis uses ordinary least squares (OLS) regression to fit:
```
delta = slope * sigma + intercept
```

Where:
- `delta` = experimental chemical shift (ppm)
- `sigma` = calculated NMR shielding (ppm)
- `slope` = scaling factor (typically -0.93 to -0.97)
- `intercept` = reference offset (typically 29-30 ppm for 1H, 171-174 ppm for 13C)

**Outlier removal:**
- Initial fit with all data points
- Identify outliers: |residual| > 3 * std(residuals)
- Refit without outliers
- Final statistics calculated from clean data

**Consistency:** 7 1H outliers and 2 13C outliers removed across all 12 solvents (identical pattern), indicating robust regression methodology.

### Factor Set Structure

Each factor entry contains:
```json
{
  "slope": -0.9340,
  "intercept": 29.79,
  "r_squared": 0.9952,
  "mae": 0.125,
  "rmsd": 0.164,
  "n_points": 335,
  "ci_slope": [-0.9410, -0.9270],
  "ci_intercept": [29.59, 29.99],
  "outliers_removed": 7
}
```

All 26 entries follow this structure - validated during merge.

### Vacuum Factors Preservation

Vacuum factors from Phase 11.2 preserved:
- 1H: slope = -0.955379905689472 (exact match)
- 13C: slope = -0.9726389627459856 (exact match)

These differ from solvent factors (more negative slope) due to lack of solvent effects on shielding.

## Next Phase Readiness

**Ready for Phase 64 - Solvent Integration:**
- All 26 factor sets available in package data
- Quality-validated for production use
- Accessible via `load_scaling_factors()` API

**Blockers:** None

**Remaining work for v2.9:**
- Update `solvent_map` in shifts.py to map user-friendly names to factor keys (Phase 64)
- Update UI to expose new solvents to users (Phase 64)
- Remove opt-in flag restrictions for new solvents (Phase 64)

## Files Modified

**Analysis code:**
- `src/qm_nmr_calc/benchmark/analysis.py` - Added 6 solvents to default list

**Package data:**
- `src/qm_nmr_calc/data/scaling_factors.json` - Merged to 26 entries (from 14)

**Benchmark output:**
- `data/benchmark/delta50/scaling_factors.json` - Generated with 24 solvent entries
- `data/benchmark/delta50/SCALING-FACTORS.md` - Updated report with all 12 solvents
- `data/benchmark/delta50/plots/*.png` - Added 24 plots for 6 new solvents

## Verification Results

**Factor count:** ✓ 26 entries in package data (13 solvents x 2 nuclei)

**R² gate:** ✓ All 12 new factors > 0.99

**1H MAE gate:** ✓ All 6 new solvents < 0.2 ppm

**13C MAE gate:** ✓ All 6 new solvents < 3.0 ppm

**Test suite:** ✓ 448 passed, 2 skipped, 0 failures (no regressions)

**Runtime load:** ✓ `load_scaling_factors()` returns 26 entries

**All success criteria met.**

## Lessons Learned

**What worked well:**
- Default solvent list approach made re-derivation simple (single-line code change)
- Merge strategy (24 new + 2 preserved) ensured consistency across all solvent factors
- Quality gates caught no issues - benchmark data quality from Phases 60-62 was excellent
- Outlier removal pattern (7 1H, 2 13C) identical across all solvents, validating robust methodology

**Observations:**
- All 6 new solvents achieved quality metrics on par with existing 6 solvents (no degradation)
- Aromatic solvents (Benzene, Toluene) continue to show best 13C performance
- Polar solvents show slightly higher 13C MAE (~2.1 ppm) vs. aromatics (~1.77 ppm), but all well below 3.0 ppm gate

**Performance:**
- Analysis runtime: ~60 seconds to derive 24 factor sets and generate plots
- Test suite runtime: ~16 minutes (no increase from baseline)
- Total execution time: 18m 38s

---

**Status:** Complete ✓
**Quality gates:** All passed ✓
**Next:** Phase 64 - Solvent Integration
