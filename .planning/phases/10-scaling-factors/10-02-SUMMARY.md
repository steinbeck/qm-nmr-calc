---
phase: 10-scaling-factors
plan: 02
subsystem: benchmark
tags: [matplotlib, report-generation, pandas, publication-quality, cli]

# Dependency graph
requires:
  - phase: 10-01-PLAN
    provides: ScalingFactor model, fit_scaling_factors(), derive_all_factors()
provides:
  - plot_regression() for scatter + residual visualization
  - plot_residual_histogram() for residual distribution
  - calculate_per_compound_stats() for per-molecule error metrics
  - generate_report() for SCALING-FACTORS.md generation
  - save_factors_json() for JSON export
  - CLI analyze command for report generation
  - Publication-quality SCALING-FACTORS.md with tables and plots
affects: [11-production-integration]

# Tech tracking
tech-stack:
  added: []  # matplotlib already in project
  patterns:
    - "Agg backend before pyplot import for headless rendering"
    - "plt.close(fig) after savefig to prevent memory leaks"
    - "orjson with OPT_INDENT_2 for human-readable JSON"

key-files:
  created:
    - data/benchmark/delta50/SCALING-FACTORS.md
    - data/benchmark/delta50/scaling_factors.json
    - data/benchmark/delta50/plots/B3LYP_CHCl3_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_CHCl3_1H_residuals.png
    - data/benchmark/delta50/plots/B3LYP_CHCl3_13C_regression.png
    - data/benchmark/delta50/plots/B3LYP_CHCl3_13C_residuals.png
    - data/benchmark/delta50/plots/B3LYP_DMSO_1H_regression.png
    - data/benchmark/delta50/plots/B3LYP_DMSO_1H_residuals.png
    - data/benchmark/delta50/plots/B3LYP_DMSO_13C_regression.png
    - data/benchmark/delta50/plots/B3LYP_DMSO_13C_residuals.png
  modified:
    - src/qm_nmr_calc/benchmark/analysis.py
    - src/qm_nmr_calc/benchmark/__main__.py
    - src/qm_nmr_calc/benchmark/__init__.py

key-decisions:
  - "2-panel regression plots: scatter + residual side-by-side"
  - "Flag compounds with mean error > 2x MAE as high-error"
  - "Sort per-compound stats by mean_error descending for easy outlier identification"

patterns-established:
  - "CLI subcommand with --factors-only flag for quick inspection"
  - "generate_report() as single entry point for all report artifacts"

# Metrics
duration: 4min
completed: 2026-01-23
---

# Phase 10 Plan 02: Report Generation Summary

**Publication-quality SCALING-FACTORS.md with regression plots achieving 0.12-2.15 ppm MAE across nuclei**

## Performance

- **Duration:** 4 min
- **Started:** 2026-01-23T13:06:33Z
- **Completed:** 2026-01-23T13:10:31Z
- **Tasks:** 3
- **Files created:** 10
- **Files modified:** 3

## Accomplishments

- plot_regression() creates 2-panel figures with scatter and residual subplots
- plot_residual_histogram() shows error distribution with MAE/RMSD annotations
- calculate_per_compound_stats() provides per-molecule error metrics for supplementary material
- generate_report() creates comprehensive SCALING-FACTORS.md (331 lines)
- save_factors_json() exports factors for external tools
- CLI analyze command with --output-dir and --factors-only options
- 8 PNG plots generated: regression + residual for each factor set

## Derived Scaling Factors

| Key | Slope | Intercept | R^2 | MAE (ppm) | n_points |
|-----|-------|-----------|-----|-----------|----------|
| B3LYP/6-311+G(2d,p)/1H/CHCl3 | -0.9375 | 29.92 | 0.9952 | 0.124 | 335 |
| B3LYP/6-311+G(2d,p)/13C/CHCl3 | -0.9497 | 172.69 | 0.9978 | 1.949 | 219 |
| B3LYP/6-311+G(2d,p)/1H/DMSO | -0.9323 | 29.73 | 0.9951 | 0.126 | 335 |
| B3LYP/6-311+G(2d,p)/13C/DMSO | -0.9429 | 171.77 | 0.9974 | 2.152 | 219 |

**MAE values compare favorably to literature benchmarks:**
- 1H: ~0.12 ppm (literature: 0.1-0.2 ppm)
- 13C: ~2.0 ppm (literature: 2-3 ppm)

## Task Commits

Each task was committed atomically:

1. **Task 1: Plot generation functions** - `db145cb` (feat)
2. **Task 2: generate_report and per-compound stats** - `a1b614b` (feat)
3. **Task 3: CLI analyze command** - `86fe6fd` (feat)
4. **Task 3: Generated report and plots** - `1224fa8` (docs)

## Files Created/Modified

**Created:**
- `data/benchmark/delta50/SCALING-FACTORS.md` - Publication-quality report (331 lines)
- `data/benchmark/delta50/scaling_factors.json` - JSON export of all factors
- `data/benchmark/delta50/plots/*.png` - 8 regression and residual plots

**Modified:**
- `src/qm_nmr_calc/benchmark/analysis.py` - Added 5 new functions
- `src/qm_nmr_calc/benchmark/__main__.py` - Added analyze subcommand
- `src/qm_nmr_calc/benchmark/__init__.py` - Exported new functions

## Decisions Made

- 2-panel regression plots show scatter and residual side-by-side for comprehensive view
- Per-compound stats sorted by mean_error descending to highlight problematic molecules
- Compounds with mean error > 2x MAE flagged in report for easy identification
- JSON export uses orjson with OPT_INDENT_2 for human readability

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - matplotlib patterns from visualization.py worked well for regression plots.

## User Setup Required

None - no external service configuration required.

## Report Highlights

The SCALING-FACTORS.md report includes:
1. **Methodology** - OLS regression with 3-sigma outlier removal
2. **Scaling Factors Table** - All 4 factor sets with 95% CIs
3. **Per-Compound Statistics** - 50 compounds per factor set for supplementary material
4. **High-Error Flagging** - Compounds with unusual errors identified
5. **Plot References** - 8 embedded PNG images
6. **Usage Examples** - Python code for applying factors

Notable high-error compounds identified:
- compound_48: 13C MAE ~21 ppm (significant outlier)
- compound_10, compound_17: Higher errors in both nuclei

## Next Phase Readiness

- Phase 10 complete: scaling factors derived and documented
- Ready for Phase 11: Production integration
- Factors ready to replace CHESHIRE TMS reference in shifts.py

---
*Phase: 10-scaling-factors*
*Completed: 2026-01-23*
