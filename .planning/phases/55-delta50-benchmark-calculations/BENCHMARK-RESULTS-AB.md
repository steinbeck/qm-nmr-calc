# DELTA50 Benchmark Results: Acetone & Benzene

**Execution Date:** 2026-02-08
**Solvents:** Acetone, Benzene
**Functional:** B3LYP / 6-311+G(2d,p)
**Molecules:** 50 (DELTA50 set from Grimblat et al. 2023)

## Summary

| Solvent | Completed | Failed | Success Rate |
|---------|-----------|--------|-------------|
| Acetone | 50/50 | 0 | 100% |
| Benzene | 50/50 | 0 | 100% |
| **Total** | **100/100** | **0** | **100%** |

## Execution Details

| Metric | Acetone | Benzene |
|--------|---------|---------|
| Total calculations | 50 | 50 |
| Average calc time | 295.3s (~4.9 min) | 295.8s (~4.9 min) |
| Min calc time | 35.1s | 35.2s |
| Max calc time | 746.8s (~12.4 min) | 753.6s (~12.6 min) |
| Total compute time | 4.2 hours | 4.2 hours |
| Processes per calc | 4 (MPI) | 4 (MPI) |

**Combined wall-clock time:** ~8.4 hours total compute (run concurrently with Plan 55-01)

## Shielding Data Quality

| Metric | Acetone | Benzene |
|--------|---------|---------|
| Total 1H data points | 342 | 342 |
| Total 13C data points | 221 | 221 |
| 1H shielding range | 21.4 to 31.6 ppm | 21.5 to 31.6 ppm |
| 13C shielding range | -56.1 to 185.7 ppm | -50.3 to 185.1 ppm |

All shielding values fall within physically expected ranges:
- 1H isotropic shielding: 21-32 ppm (typical for organic H atoms)
- 13C isotropic shielding: -56 to 186 ppm (negative values for carbonyl/unsaturated C, positive for aliphatic C)

Acetone and Benzene produce very similar shielding values, as expected -- solvent effects are real but moderate for COSMO implicit solvation. Benzene shows a slightly narrower 13C range (-50.3 vs -56.1 ppm minimum), consistent with benzene's lower dielectric constant providing less electrostatic stabilization of polar functional groups.

## Failures

None -- all 100 calculations completed successfully.

## Bug Fixes Applied

Two bug fixes were required during the pilot run phase (commit df3ca05):

1. **NWChem scratch file isolation** -- Sequential benchmark calculations wrote `molecule.movecs` and other scratch files to the shared project root directory. When the next calculation started, stale scratch files caused "could not read mo vectors" errors. Fixed by setting `cwd` to the per-calculation scratch directory in `run_nwchem()`.

2. **CPHF convergence with COSMO** -- NMR property calculations failed with `cphf_solve2: SCF residual greater than 1d-2`. Fixed by adding the `direct` keyword to the DFT block in the shielding input generator, forcing on-the-fly integral evaluation required for stable CPHF convergence with COSMO solvation.

Both fixes also benefited the concurrent Plan 55-01 (Methanol & Water) calculations.

## Phase 55 Success Criteria Verification

With both Plan 55-01 and 55-02 complete, all Phase 55 success criteria are satisfied:

| # | Criterion | Status |
|---|-----------|--------|
| 1 | Methanol benchmark: 50 completed NWChem outputs with NMR shielding tensors | 50/50 |
| 2 | Water benchmark: 50 completed NWChem outputs with NMR shielding tensors | 50/50 |
| 3 | Acetone benchmark: 50 completed NWChem outputs with NMR shielding tensors | 50/50 |
| 4 | Benzene benchmark: 50 completed NWChem outputs with NMR shielding tensors | 50/50 |
| 5 | No COSMO convergence errors across all 200 calculations | 0 failures |

**Total: 200/200 calculations successful across all 4 solvents.**

## Data Location

Results: `data/benchmark/results/compound_XX/B3LYP_{Acetone|Benzene}/shifts.json`

Each shifts.json contains:
- `shielding_data`: Raw isotropic shielding tensors (used for scaling factor derivation)
  - `index`: Atom indices (1-based)
  - `atom`: Element symbols
  - `shielding`: Isotropic shielding values in ppm
- `h1_shifts` / `c13_shifts`: Empty arrays (no scaling factors yet for these solvents)

## Next Steps

- Phase 56 derives scaling factors from all 4 solvents' shielding data via linear regression
