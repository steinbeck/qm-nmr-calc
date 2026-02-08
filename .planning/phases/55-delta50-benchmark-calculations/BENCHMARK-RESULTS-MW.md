# DELTA50 Benchmark Results: Methanol & Water

**Execution Date:** 2026-02-08
**Solvents:** Methanol, Water
**Functional:** B3LYP / 6-311+G(2d,p)
**Molecules:** 50 (DELTA50 set from Grimblat et al. 2023)

## Summary

| Solvent | Completed | Failed | Success Rate |
|---------|-----------|--------|-------------|
| Methanol | 50/50 | 0 | 100% |
| Water | 50/50 | 0 | 100% |
| **Total** | **100/100** | **0** | **100%** |

## Execution Details

| Metric | Methanol | Water |
|--------|----------|-------|
| Total calculations | 50 | 50 |
| Average calc time | 275.1s (~4.6 min) | 275.5s (~4.6 min) |
| Min calc time | 20.6s | 19.7s |
| Max calc time | 735.9s (~12.3 min) | 732.2s (~12.2 min) |
| Total compute time | 4.2 hours | 4.2 hours |
| Processes per calc | 4 (MPI) | 4 (MPI) |

**Combined wall-clock time:** ~8.4 hours total compute (run concurrently)

## Shielding Data Quality

| Metric | Methanol | Water |
|--------|----------|-------|
| Total 1H data points | 342 | 342 |
| Total 13C data points | 221 | 221 |
| 1H shielding range | 21.4 to 31.6 ppm | 21.4 to 31.6 ppm |
| 13C shielding range | -56.4 to 185.7 ppm | -56.7 to 185.7 ppm |

All shielding values fall within physically expected ranges:
- 1H isotropic shielding: 21-32 ppm (typical for organic H atoms)
- 13C isotropic shielding: -57 to 186 ppm (negative values for carbonyl/unsaturated C, positive for aliphatic C)

Methanol and Water produce very similar shielding values, as expected -- solvent effects are real but moderate for COSMO implicit solvation.

## Failures

None -- all 100 calculations completed successfully.

## Bug Fix Applied

During the pilot run, all calculations failed with `cphf_solve2: SCF residual greater than 1d-2` (CPHF convergence error). This was fixed by adding the `direct` keyword to the DFT block in the NMR shielding input generator (`input_gen.py`), which forces on-the-fly integral evaluation required for stable CPHF convergence with COSMO solvation. Commit: df3ca05.

## Data Location

Results: `data/benchmark/results/compound_XX/B3LYP_{Methanol|Water}/shifts.json`

Each shifts.json contains:
- `shielding_data`: Raw isotropic shielding tensors (used for scaling factor derivation)
  - `index`: Atom indices (1-based)
  - `atom`: Element symbols
  - `shielding`: Isotropic shielding values in ppm
- `h1_shifts` / `c13_shifts`: Empty arrays (no scaling factors yet for these solvents)

## Next Steps

- Plan 55-02 runs Acetone and Benzene benchmarks (already completed concurrently)
- Phase 56 derives scaling factors from all 4 solvents' shielding data
