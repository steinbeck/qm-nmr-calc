---
phase: 56-scaling-factor-derivation
verified: 2026-02-09T08:57:08Z
status: passed
score: 5/5 must-haves verified
---

# Phase 56: Scaling Factor Derivation Verification Report

**Phase Goal:** OLS-derived scaling factors for all 4 new solvents pass quality gates and are stored in package data

**Verified:** 2026-02-09T08:57:08Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | analyze command produces scaling factors for all 6 B3LYP solvents (CHCl3, DMSO, Methanol, Water, Acetone, Benzene) | ✓ VERIFIED | analysis.py line 213 contains all 6 solvents in default list; delta50/scaling_factors.json contains 12 entries (6 solvents x 2 nuclei) |
| 2 | All 8 new factor sets (4 solvents x 2 nuclei) have R-squared > 0.99 | ✓ VERIFIED | Methanol: 1H=0.9951, 13C=0.9974; Water: 1H=0.9951, 13C=0.9974; Acetone: 1H=0.9951, 13C=0.9975; Benzene: 1H=0.9950, 13C=0.9980 |
| 3 | 1H MAE is below 0.2 ppm for each of the 4 new solvents | ✓ VERIFIED | Methanol: 0.126 ppm, Water: 0.127 ppm, Acetone: 0.126 ppm, Benzene: 0.128 ppm — all < 0.2 ppm |
| 4 | 13C MAE is below 3.0 ppm for each of the 4 new solvents | ✓ VERIFIED | Methanol: 2.139 ppm, Water: 2.161 ppm, Acetone: 2.117 ppm, Benzene: 1.761 ppm — all < 3.0 ppm |
| 5 | Package scaling_factors.json contains 14 total factor sets (6 existing + 8 new) | ✓ VERIFIED | src/qm_nmr_calc/data/scaling_factors.json has 14 entries: 12 solvent factors (6 solvents x 2 nuclei) + 2 vacuum factors |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/qm_nmr_calc/benchmark/analysis.py` | Updated default solvents list | ✓ VERIFIED | Line 213: `solvents = ["CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene"]` — substantive (429 lines), wired (imported by __main__.py, aggregate_regression_data reads from results_dir) |
| `src/qm_nmr_calc/data/scaling_factors.json` | Package scaling factors with all 14 factor sets | ✓ VERIFIED | 14 entries, all with slope/intercept/r_squared/mae/rmsd/n_points/ci_slope/ci_intercept/outliers_removed — substantive (240 lines), wired (loaded by shifts.py get_scaling_factor) |
| `data/benchmark/delta50/scaling_factors.json` | Analysis output with all 12 derived factor sets | ✓ VERIFIED | 12 entries for 6 solvents x 2 nuclei — substantive (206 lines), wired (copied to package data, all entries match) |
| `data/benchmark/delta50/SCALING-FACTORS.md` | Updated report with all solvents | ✓ VERIFIED | Contains tables for all 6 solvents (Methanol, Water, Acetone, Benzene, CHCl3, DMSO) — substantive (714 lines), includes quality metrics and per-compound statistics |

### Key Link Verification

| From | To | Via | Status | Details |
|------|--|----|--------|---------|
| analysis.py | benchmark results | aggregate_regression_data | ✓ WIRED | Line 67: `shifts_file = results_dir / mol_id / f"{functional}_{solvent}" / "shifts.json"` — reads shielding values from results directory |
| delta50/scaling_factors.json | package scaling_factors.json | copy/merge | ✓ WIRED | All 12 delta50 entries present in package data with identical values; 2 vacuum entries added separately |
| shifts.py | package scaling_factors.json | load_scaling_factors | ✓ WIRED | Runtime test successful: get_scaling_factor returns valid factors for all 4 new solvents; 14 total factors loaded |

### Requirements Coverage

No specific requirements mapped to this phase in REQUIREMENTS.md. Phase delivers enabling infrastructure for NMR prediction accuracy.

### Anti-Patterns Found

None. Clean implementation with no TODO/FIXME comments, placeholders, or stub patterns detected.

### Quality Metrics Detail

**New Solvents (8 factor sets):**

| Solvent | Nucleus | R² | MAE (ppm) | RMSD (ppm) | n_points | Quality Gates |
|---------|---------|-----|-----------|------------|----------|---------------|
| Methanol | 1H | 0.9951 | 0.126 | 0.166 | 335 | ✓ R²>0.99, MAE<0.2 |
| Methanol | 13C | 0.9974 | 2.139 | 2.901 | 219 | ✓ R²>0.99, MAE<3.0 |
| Water | 1H | 0.9951 | 0.127 | 0.166 | 335 | ✓ R²>0.99, MAE<0.2 |
| Water | 13C | 0.9974 | 2.161 | 2.927 | 219 | ✓ R²>0.99, MAE<3.0 |
| Acetone | 1H | 0.9951 | 0.126 | 0.165 | 335 | ✓ R²>0.99, MAE<0.2 |
| Acetone | 13C | 0.9975 | 2.117 | 2.875 | 219 | ✓ R²>0.99, MAE<3.0 |
| Benzene | 1H | 0.9950 | 0.128 | 0.167 | 335 | ✓ R²>0.99, MAE<0.2 |
| Benzene | 13C | 0.9980 | 1.761 | 2.538 | 219 | ✓ R²>0.99, MAE<3.0 |

**All 8 new factor sets pass quality gates.**

**Existing Solvents (re-derived):**

| Solvent | Nucleus | R² | MAE (ppm) | n_points |
|---------|---------|-----|-----------|----------|
| CHCl3 | 1H | 0.9952 | 0.124 | 335 |
| CHCl3 | 13C | 0.9978 | 1.949 | 219 |
| DMSO | 1H | 0.9951 | 0.126 | 335 |
| DMSO | 13C | 0.9974 | 2.152 | 219 |

**Vacuum (preserved from Phase 11.2):**

| Nucleus | Slope | Intercept | R² | MAE (ppm) | n_points |
|---------|-------|-----------|-----|-----------|----------|
| 1H | -0.9554 | 30.5446 | 0.9934 | 0.148 | 336 |
| 13C | -0.9726 | 175.7109 | 0.9980 | 1.739 | 219 |

### Diagnostic Outputs

**Plots generated:** 28 total
- Regression scatter plots: 14 (7 contexts x 2 nuclei)
- Residual histograms: 14 (7 contexts x 2 nuclei)
- Contexts: CHCl3, DMSO, Methanol, Water, Acetone, Benzene, vacuum

All plots exist in `data/benchmark/delta50/plots/` with reasonable file sizes (40-160 KB).

**Report generated:**
- `data/benchmark/delta50/SCALING-FACTORS.md`: 714 lines
- Contains methodology, factor tables, statistical summaries, per-compound error analysis

### Runtime Validation

**Scaling factor loading test:**
```python
from qm_nmr_calc.shifts import get_scaling_factor
get_scaling_factor('B3LYP', '6-311+G(2d,p)', '1H', 'Methanol')
# Returns: {'slope': -0.9326, 'intercept': 29.7450, 'r_squared': 0.9951, 'mae': 0.126, ...}
```

✓ All 4 new solvents accessible via shifts.py API
✓ load_scaling_factors() returns 14 entries
✓ Vacuum factors preserved exactly

### Verification Summary

**Phase goal achieved:** OLS-derived scaling factors for all 4 new solvents pass quality gates and are stored in package data.

**Evidence:**
1. All 8 new factor sets derived with R² > 0.99
2. All 1H MAE < 0.2 ppm, all 13C MAE < 3.0 ppm
3. Package data contains 14 total factor sets
4. Delta50 analysis output contains 12 solvent factor sets
5. Runtime loading successful for all new solvents
6. Comprehensive regression report and plots generated
7. Default solvents list updated to include all 6 solvents

**No gaps found.** All must-haves verified. Phase ready for downstream integration.

---

*Verified: 2026-02-09T08:57:08Z*
*Verifier: Claude (gsd-verifier)*
