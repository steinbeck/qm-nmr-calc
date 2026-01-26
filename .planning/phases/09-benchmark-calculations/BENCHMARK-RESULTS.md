# DELTA50 Benchmark Results

**Execution Date:** 2026-01-22
**Total Duration:** ~2.2 hours (main run) + ~15 min (compound_14 re-run)

## Summary

| Metric | Value |
|--------|-------|
| Total Molecules | 50 |
| Functionals Run | B3LYP |
| Solvents | CHCl3, DMSO |
| Total Calculations | 100 |
| Completed | 100 |
| Failed | 0 |
| Success Rate | 100% |

**Note:** Only B3LYP functional was included in this benchmark run. WP04 was not executed as originally planned in the 200-calculation scope.

## Results by Method

| Functional | Solvent | Completed | Failed |
|------------|---------|-----------|--------|
| B3LYP | CHCl3 | 50 | 0 |
| B3LYP | DMSO | 50 | 0 |
| WP04 | CHCl3 | 0 | - |
| WP04 | DMSO | 0 | - |

## Shift Statistics

### 1H Chemical Shifts
| Statistic | Value (ppm) |
|-----------|-------------|
| Count | 684 |
| Minimum | -1.50 |
| Maximum | 9.27 |
| Average | 1.91 |

### 13C Chemical Shifts
| Statistic | Value (ppm) |
|-----------|-------------|
| Count | 442 |
| Minimum | -10.72 |
| Maximum | 239.11 |
| Average | 84.46 |

The shift ranges are within expected bounds for typical organic molecules:
- 1H: -2 to 15 ppm typical range
- 13C: -10 to 250 ppm typical range

## Issues Encountered and Resolved

### compound_14 (2-butyne) - AUTOZ Failure

**Initial Error:** NWChem's AUTOZ (automatic internal coordinate generation) failed for this linear alkyne molecule:
```
!! There are insufficient internal variables: expected 23 got 24
!! Either AUTOZ failed or your geometry has changed so much that the
!! coordinates should be regenerated.
geom_binvr: #indep variables incorrect
```

**Root Cause:** Linear molecules (like alkynes) can cause AUTOZ to fail when it cannot generate proper internal coordinates.

**Fix Applied:** Added optional `noautoz` parameter to `generate_optimization_input()` and `generate_shielding_input()` functions in `src/qm_nmr_calc/nwchem/input_gen.py`. When enabled, this forces NWChem to use Cartesian coordinates instead of internal coordinates.

**Important:** `noautoz` should NOT be used routinely - only as a fallback for molecules that fail with AUTOZ errors. The default remains `noautoz=False`.

**Resolution:** compound_14 was successfully re-run with `noautoz=True`, completing both CHCl3 and DMSO calculations.

## Data Location

Results stored in: `data/benchmark/results/`

### Directory Structure
```
data/benchmark/results/
├── compound_01/
│   ├── B3LYP_CHCl3/
│   │   └── shifts.json
│   └── B3LYP_DMSO/
│       └── shifts.json
├── compound_02/
│   └── ...
└── ...
```

### File Format
Each `shifts.json` contains:
```json
{
    "molecule_id": "compound_XX",
    "functional": "B3LYP",
    "solvent": "CHCl3|DMSO",
    "h1_shifts": [1.23, 2.34, ...],
    "c13_shifts": [45.6, 78.9, ...],
    "shielding_data": {...}
}
```

## Benchmark Runner Status

The `status.json` file shows historical failure records from before the compound_14 fix was applied. The actual data is complete (100/100 shifts.json files present).

## Performance Metrics

| Metric | Value |
|--------|-------|
| Average calculation time | ~109 seconds |
| Total run time | ~2.2 hours |
| Processes used | 8 (parallel) |

## Next Steps

Phase 10 will use these results to:
1. Load calculated shifts from `shifts.json` files
2. Load experimental shifts from DELTA50 dataset
3. Derive scaling factors via linear regression
4. Evaluate prediction accuracy (MAE, RMSE, R²)

### Note on WP04 Functional

The WP04 functional calculations were not included in this run. If WP04 results are needed for comparison:
1. Update the benchmark runner to queue WP04 calculations
2. Run: `python -m qm_nmr_calc.benchmark run --headless --processes 8`
3. This would add another ~100 calculations (~2 hours)
