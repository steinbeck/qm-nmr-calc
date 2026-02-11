# DELTA50 Benchmark Results: Pyridine & THF

**Execution Date:** 2026-02-10 to 2026-02-11
**Solvents:** Pyridine, THF
**Functional:** B3LYP / 6-311+G(2d,p)
**Molecules:** 50 (DELTA50 set)

## Summary

| Solvent | Completed | Failed | Success Rate |
|---------|-----------|--------|-------------|
| Pyridine | 50/50 | 0 | 100% |
| THF | 50/50 | 0 | 100% |
| **Total** | **100/100** | **0** | **100%** |

## Timing

- Total runtime: ~10.5 hours (pilot + full headless run)
- Pilot run: 10 calculations (5 Pyridine, 5 THF) - all successful
- Full run: 90 calculations (45 Pyridine, 45 THF) - all successful
- COSMO parameters:
  - Pyridine (dielectric=12.978)
  - THF (dielectric=7.4257)

## Failures

None. All 100 calculations completed successfully with valid shielding data.

## Data Quality Verification

Spot-checked results for compounds 01, 25, and 50:
- All contain valid `shielding_data` arrays
- Proper atom type identification (H, C, N, O)
- Shielding tensor values in expected ranges
- Empty `h1_shifts` and `c13_shifts` arrays (awaiting scaling factors)

## Data Location

Results: `data/benchmark/results/compound_XX/B3LYP_{Pyridine|THF}/shifts.json`

Each shifts.json contains:
- `shielding_data`: Raw isotropic shielding tensors (used for scaling factor derivation)
- `h1_shifts` / `c13_shifts`: Empty (no scaling factors yet for these solvents)

## Requirements Satisfied

- **BENCH-02**: All 50 Pyridine DELTA50 molecules calculated successfully with COSMO solvation ✓
- **BENCH-03**: All 50 THF DELTA50 molecules calculated successfully with COSMO solvation ✓
- Calculated shielding data extracted and stored for both solvents ✓

## Next Steps

- Phase 61 runs Toluene and DCM benchmarks (2 more solvents, 100 calculations)
- Phase 62 runs DMF and Acetone benchmarks (final 2 solvents, 100 calculations)
- Phase 63 derives scaling factors from all 6 solvents' combined shielding data
