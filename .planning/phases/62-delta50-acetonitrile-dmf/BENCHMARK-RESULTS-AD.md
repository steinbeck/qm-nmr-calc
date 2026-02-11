# DELTA50 Benchmark Results: Acetonitrile & DMF

**Execution Date:** 2026-02-10 to 2026-02-11
**Solvents:** Acetonitrile, DMF
**Functional:** B3LYP / 6-311+G(2d,p)
**Molecules:** 50 (DELTA50 set)

## Summary

| Solvent | Completed | Failed | Success Rate |
|---------|-----------|--------|-------------|
| Acetonitrile | 50/50 | 0 | 100% |
| DMF | 50/50 | 0 | 100% |
| **Total** | **100/100** | **0** | **100%** |

## Timing

- Total runtime: ~10.5 hours (shared with Phase 60 - all 6 solvents run in parallel)
- Calculations run in parallel with Phase 60 using --processes 30
- Total compute time shared with Phase 60: ~10.5 hours for all 6 solvents (Pyridine, THF, Toluene, DCM, Acetonitrile, DMF)
- COSMO parameters:
  - Acetonitrile (dielectric~37.5, NWChem COSMO name "acetntrl")
  - DMF (dielectric~36.7)

## Failures

None. All 100 calculations completed successfully with valid shielding data.

## Data Quality Verification

Spot-checked results for compounds 01, 25, and 50:
- compound_01/Acetonitrile: 3H, 1C, 7 total atoms
- compound_01/DMF: 3H, 1C, 7 total atoms
- compound_25/Acetonitrile: 4H, 4C, 10 total atoms
- compound_25/DMF: 4H, 4C, 10 total atoms
- compound_50/Acetonitrile: 6H, 4C, 11 total atoms
- compound_50/DMF: 6H, 4C, 11 total atoms
- All contain valid `shielding_data` arrays
- Proper atom type identification (H, C, N, O)
- Shielding tensor values in expected ranges
- Empty `h1_shifts` and `c13_shifts` arrays (awaiting scaling factors)

## Data Location

Results: `data/benchmark/results/compound_XX/B3LYP_{Acetonitrile|DMF}/shifts.json`

Each shifts.json contains:
- `shielding_data`: Raw isotropic shielding tensors (used for scaling factor derivation)
- `h1_shifts` / `c13_shifts`: Empty (no scaling factors yet for these solvents)

## Requirements Satisfied

- **BENCH-06**: All 50 acetonitrile DELTA50 molecules calculated successfully with COSMO solvation (using acetntrl COSMO name) ✓
- **BENCH-07**: All 50 DMF DELTA50 molecules calculated successfully with COSMO solvation ✓
- Calculated shielding data extracted and stored for both solvents ✓

## Next Steps

- Phase 63 derives scaling factors from all 6 solvents' combined shielding data (300 total calculations: Pyridine, THF, Toluene, DCM, Acetonitrile, DMF)
