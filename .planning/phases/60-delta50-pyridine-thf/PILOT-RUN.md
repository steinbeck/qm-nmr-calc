# Pilot Run - Pyridine and THF Benchmark

**Date:** 2026-02-10
**Phase:** 60-delta50-pyridine-thf
**Plan:** 60-01

## Execution Summary

**Molecules tested:** compound_01, compound_02, compound_03, compound_04, compound_05
**Solvents:** Pyridine, THF
**Total calculations:** 10 (5 molecules x 2 solvents)

**Status:** All 10 calculations completed successfully
**Failures:** 0
**Average calculation time:** 38.4 seconds

## Results Validation

All pilot results verified:

- **compound_01/Pyridine:** H=3 (σ=27.1-27.4 ppm), C=1 (σ=115.1 ppm)
- **compound_01/THF:** H=3 (σ=27.2-27.4 ppm), C=1 (σ=115.3 ppm)
- **compound_03/Pyridine:** H=4 (σ=21.4-29.9 ppm), C=2 (σ=-32.9-147.6 ppm)
- **compound_03/THF:** H=4 (σ=21.4-29.9 ppm), C=2 (σ=-32.1-147.7 ppm)
- **compound_05/Pyridine:** H=3 (σ=29.7 ppm), C=2 (σ=54.7-182.2 ppm)
- **compound_05/THF:** H=3 (σ=29.7 ppm), C=2 (σ=55.4-182.2 ppm)

## Quality Checks

- ✓ All shifts.json files contain `shielding_data` with atom entries
- ✓ H shielding values in physically reasonable ranges (~21-30 ppm)
- ✓ C shielding values in physically reasonable ranges (~-33 to 182 ppm)
- ✓ Empty h1_shifts/c13_shifts arrays (expected - no scaling factors yet)
- ✓ No COSMO convergence errors
- ✓ No SCF convergence failures

## COSMO Solvation Validation

Both Pyridine and THF solvents working correctly with NWChem COSMO:
- Pyridine mapped to "pyridine" COSMO solvent
- THF mapped to "thf" COSMO solvent
- All calculations converged without solvation errors

## Next Steps

Ready to proceed with full headless run for remaining 45 molecules (90 calculations).
Estimated completion time: ~8-10 hours based on 38.4s average.
