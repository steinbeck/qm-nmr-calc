# DELTA50 Benchmark Results: Toluene & DCM

**Execution Date:** 2026-02-10 to 2026-02-11
**Solvents:** Toluene, DCM
**Functional:** B3LYP / 6-311+G(2d,p)
**Molecules:** 50 (DELTA50 set)

## Summary

| Solvent | Completed | Failed | Success Rate |
|---------|-----------|--------|-------------|
| Toluene | 50/50 | 0 | 100% |
| DCM | 50/50 | 0 | 100% |
| **Total** | **100/100** | **0** | **100%** |

## Timing

- Calculations run in parallel with Phase 60 using --processes 30
- Runtime: ~10.5 hours (concurrent with Pyridine/THF benchmark)
- COSMO parameters:
  - Toluene (dielectric=2.3741)
  - DCM (dielectric=8.930)

## Failures

None. All 100 calculations completed successfully with valid shielding data.

## Data Quality Verification

Spot-checked results for compounds 01, 25, and 50:
- compound_01/Toluene: 3H, 1C, 7 total atoms
- compound_01/DCM: 3H, 1C, 7 total atoms
- compound_25/Toluene: 4H, 4C, 10 total atoms
- compound_25/DCM: 4H, 4C, 10 total atoms
- compound_50/Toluene: 6H, 4C, 11 total atoms
- compound_50/DCM: 6H, 4C, 11 total atoms

All contain valid `shielding_data` arrays with proper atom type identification (H, C, N, O) and shielding tensor values in expected ranges. Empty `h1_shifts` and `c13_shifts` arrays (awaiting scaling factors).

## Data Location

Results: `data/benchmark/results/compound_XX/B3LYP_{Toluene|DCM}/shifts.json`

Each shifts.json contains:
- `shielding_data`: Raw isotropic shielding tensors (used for scaling factor derivation)
- `h1_shifts` / `c13_shifts`: Empty (no scaling factors yet for these solvents)

## Requirements Satisfied

- **BENCH-04**: All 50 Toluene DELTA50 molecules calculated successfully with COSMO solvation ✓
- **BENCH-05**: All 50 DCM DELTA50 molecules calculated successfully with COSMO solvation ✓
- Calculated shielding data extracted and stored for both solvents ✓

## Next Steps

- Phase 62 runs DMF and Acetone benchmarks (final 2 solvents, 100 calculations)
- Phase 63 derives scaling factors from all 6 solvents' combined shielding data (300 total calculations: Pyridine, THF, Toluene, DCM, DMF, Acetone)
