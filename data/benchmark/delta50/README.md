# DELTA50 Benchmark Dataset

## Citation

Grimblat, N.; Puente, Á.R.; Sarotti, A.M. "Computational Chemistry Driven Solution to Rubrofusarin Gentiobioside Structure". *Molecules* **2023**, *28*(6), 2449.

**DOI:** [10.3390/molecules28062449](https://doi.org/10.3390/molecules28062449)

## Description

The DELTA50 test set consists of 50 organic molecules with known experimental 1H and 13C NMR chemical shifts. This benchmark dataset is used to validate the accuracy of computational NMR prediction methods.

## Data Sources

### Molecule Structures (XYZ Files)

**Required:** 50 XYZ coordinate files in `molecules/` directory

**Source:** Supporting Information from the publication above
- Location: https://www.mdpi.com/1420-3049/28/6/2449/s1
- File format: XYZ coordinates (Cartesian)
- Naming: `molecule_01.xyz` through `molecule_50.xyz`

### Experimental Chemical Shifts

**Required:** `experimental_shifts.json`

**Source:** Excel file in Supporting Information containing experimental 1H and 13C chemical shifts
- Conditions: CDCl3, 298K (standard for DELTA50)
- Format: Parsed into JSON with atom assignments

## Usage

To use this benchmark data in Python:

```python
from qm_nmr_calc.benchmark import load_delta50_molecules, load_experimental_shifts

# Load molecular structures
molecules = load_delta50_molecules()  # Dict[str, Tuple[Mol, str]]

# Load experimental shifts
shifts = load_experimental_shifts()  # ExperimentalShifts object
```

## Manual Setup Required

**Important:** Due to download restrictions on MDPI's SI portal, the data files must be obtained manually:

1. Visit: https://www.mdpi.com/1420-3049/28/6/2449/s1
2. Download the Supporting Information ZIP file
3. Extract XYZ coordinate files to `molecules/` directory
4. Parse experimental shifts from Excel file to `experimental_shifts.json`

The molecular structures and experimental data are published under the Creative Commons Attribution (CC BY) license as part of the MDPI Molecules journal.
