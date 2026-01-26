# DELTA50 Benchmark Dataset

## Citation

**DELTA50 Dataset Paper:** Grimblat, N.; Gavin, J. A.; Hernandez Daranas, A.; Sarotti, A. M. *Combining the Power of J Coupling and DP4 Analysis on Stereochemical Assignments: The J-DP4 Method.* Molecules **2023**, 28(6), 2449.

**DOI:** [10.3390/molecules28062449](https://doi.org/10.3390/molecules28062449)

**Original DP4 Paper:** Grimblat, N.; Zanardi, M. M.; Sarotti, A. M. *Beyond DP4: An Improved Probability for the Stereochemical Assignment of Isomeric Compounds Using Quantum Chemical Calculations of NMR Shifts.* J. Org. Chem. **2015**, 80, 12526-12534.

## Description

The DELTA50 dataset consists of 50 small organic molecules used as a benchmark for NMR chemical shift calculations. The dataset includes:

- 50 molecules with diverse functional groups
- B3LYP/6-31G(d) optimized geometries (in vacuo)
- Experimental 1H and 13C chemical shifts measured in CDCl3 at 298K
- Chemical shifts referenced to TMS at 0.00 ppm

## Contents

### Molecules (`molecules/`)

50 XYZ structure files named `compound_01.xyz` through `compound_50.xyz`:

| # | Name | # | Name |
|---|------|---|------|
| 1 | Nitromethane | 26 | Cyclopentane |
| 2 | Nitroethane | 27 | 2-Methyl-2-butene |
| 3 | Acetaldehyde | 28 | Methyl t-butyl ether (MTBE) |
| 4 | Oxirane | 29 | Tetrahydropyran (THP) |
| 5 | Acetonitrile | 30 | N-Methylpyrrolidine |
| 6 | Cyclopropane | 31 | Cyclopentanone |
| 7 | Acetone | 32 | Pivalonitrile |
| 8 | Oxetane | 33 | Cyclopent-2-en-1-one |
| 9 | Methyl acetate | 34 | N-Methylpyrrole |
| 10 | N,N-Dimethylformamide (DMF) | 35 | Pyridine |
| 11 | Propionitrile | 36 | t-Butylethylene |
| 12 | Isoxazole | 37 | Cyclohexane |
| 13 | Isobutylene | 38 | t-Butylacetylene |
| 14 | 2-Butyne | 39 | Benzene |
| 15 | t-Butyl nitrate | 40 | N-Methylpiperidine |
| 16 | Tetrahydrofuran (THF) | 41 | Cyclohexanone |
| 17 | N,N-Dimethylacetamide (DMAc) | 42 | Cyclohex-2-en-1-one |
| 18 | Cyclobutanone | 43 | Fluorobenzene |
| 19 | Butyrolactone | 44 | Nitrobenzene |
| 20 | Isobutyronitrile | 45 | 1,4-Benzoquinone |
| 21 | Furan | 46 | Toluene |
| 22 | 3-Butyn-2-one | 47 | Norbornadiene |
| 23 | Pyrimidine | 48 | Anisole |
| 24 | 1,4-Pyrazine | 49 | Maleic anhydride |
| 25 | Pyridazine | 50 | 2,5-Dihydrofuran |

### Experimental Shifts (`experimental_shifts.json`)

JSON file containing:
- Source citation and experimental conditions
- Per-molecule data:
  - Compound name
  - Unique 1H and 13C chemical shifts (ppm)
  - Atom-index-to-shift mappings for both nuclei

## Usage

```python
from qm_nmr_calc.benchmark import load_delta50_molecules, load_experimental_shifts

# Load all molecule structures
molecules = load_delta50_molecules()  # Dict[str, (RDKit Mol, XYZ block)]

# Load experimental shifts
shifts = load_experimental_shifts()  # ExperimentalShifts model
print(f"Loaded {len(shifts.molecules)} molecules")

# Access data for a specific molecule
mol_data = shifts.molecules["compound_39"]  # Benzene
print(f"Benzene 13C shifts: {mol_data.c13_shifts}")  # [128.34]
print(f"Benzene 1H shifts: {mol_data.h1_shifts}")    # [7.36]
```

## Notes

- Geometries are from the original publication's Supporting Information
- The dominant conformer (>98% Boltzmann population) is provided for each compound
- Heteroatom (N, O, F) experimental shifts are not available (marked as "-" in original data)
- All coordinates are in Angstroms

## License

Data extracted from published Supporting Information for research purposes.
