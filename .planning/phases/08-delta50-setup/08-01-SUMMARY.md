---
phase: 08-delta50-setup
plan: 01
subsystem: benchmark
tags: [delta50, nmr, benchmark, data]

dependency-graph:
  requires: []
  provides: [delta50-structures, experimental-shifts, benchmark-loader]
  affects: [08-02, 08-03, 09-benchmark-execution]

tech-stack:
  added: []
  patterns: [pdf-parsing, xyz-format, json-data-storage]

key-files:
  created:
    - data/benchmark/delta50/molecules/compound_01.xyz through compound_50.xyz
    - data/benchmark/delta50/experimental_shifts.json
    - data/benchmark/delta50/README.md
    - src/qm_nmr_calc/benchmark/__init__.py
    - src/qm_nmr_calc/benchmark/models.py
    - src/qm_nmr_calc/benchmark/data_loader.py
    - scripts/generate_delta50_data.py
  modified: []

decisions:
  - id: compound-naming
    choice: "compound_XX.xyz naming convention"
    reason: "Clearer than molecule_XX, matches original compound numbering from paper"
  - id: json-shift-structure
    choice: "Both unique lists and atom-index mappings stored"
    reason: "Supports both quick lookup and atom-level analysis"
  - id: pdf-extraction
    choice: "Manual parsing from PDF to Python dict"
    reason: "SI only available as PDF, not machine-readable files"

metrics:
  duration: "9 min"
  completed: "2026-01-21"
---

# Phase 08 Plan 01: DELTA50 Benchmark Data Acquisition Summary

**One-liner:** Extract 50 molecule structures and experimental shifts from DELTA50 paper SI PDF into XYZ files and JSON format with Python loader functions.

## What Was Built

### 1. XYZ Structure Files (50 files)
- Created `data/benchmark/delta50/molecules/compound_01.xyz` through `compound_50.xyz`
- B3LYP/6-31G(d) optimized geometries from the publication SI
- Standard XYZ format with atom count, title line, and Cartesian coordinates
- Molecules range from 6 atoms (acetonitrile) to 20 atoms (N-methylpiperidine)

### 2. Experimental Shifts JSON
- Created `data/benchmark/delta50/experimental_shifts.json`
- Contains 1H and 13C chemical shifts for all 50 molecules
- Includes atom-index-to-shift mappings for detailed analysis
- Source citation and experimental conditions documented

### 3. Benchmark Module
- `src/qm_nmr_calc/benchmark/models.py`: Pydantic models (MoleculeData, ExperimentalShifts, BenchmarkResult)
- `src/qm_nmr_calc/benchmark/data_loader.py`: Functions to load data (load_delta50_molecules, load_experimental_shifts)
- `src/qm_nmr_calc/benchmark/__init__.py`: Public exports

### 4. Generation Script
- `scripts/generate_delta50_data.py`: Reproducible data generation from parsed PDF content
- Contains all 50 compounds with coordinates and shifts as Python dicts
- Can regenerate XYZ files and JSON if needed

## Key Implementation Details

### Data Source
The DELTA50 data was extracted from the paper's Supporting Information PDF, which contained:
- Tables with atom indices, element types, x/y/z coordinates, and experimental shifts
- B3LYP/6-31G(d) optimized geometries in vacuo
- Experimental shifts measured in CDCl3 at 298K, referenced to TMS at 0.00 ppm

### Shift Organization
```python
# Example entry in experimental_shifts.json
"compound_39": {
    "name": "Benzene",
    "compound_number": 39,
    "h1_shifts": [7.36],           # Unique values
    "c13_shifts": [128.34],        # Unique values
    "h1_assignments": {"7": 7.36, "8": 7.36, ...},  # Atom index -> shift
    "c13_assignments": {"1": 128.34, "2": 128.34, ...},
    "num_h_atoms": 6,
    "num_c_atoms": 6
}
```

### Usage Example
```python
from qm_nmr_calc.benchmark import load_delta50_molecules, load_experimental_shifts

# Load structures
molecules = load_delta50_molecules()  # {id: (RDKit Mol, XYZ block)}
print(f"Loaded {len(molecules)} molecules")

# Load experimental data
shifts = load_experimental_shifts()
benzene = shifts.molecules["compound_39"]
print(f"Benzene 13C: {benzene.c13_shifts} ppm")
```

## Deviations from Plan

None - plan executed exactly as written. The plan was updated to use PDF extraction since the SI data was embedded in a PDF rather than available as downloadable files.

## Commits

| Hash | Message |
|------|---------|
| 0d6f928 | feat(08-01): add 50 DELTA50 molecule XYZ structures |
| a715f83 | feat(08-01): add DELTA50 experimental shifts JSON |
| 0083be4 | feat(08-01): create benchmark module with data loader |

## Verification Results

```
$ ls data/benchmark/delta50/molecules/*.xyz | wc -l
50

$ python -c "from qm_nmr_calc.benchmark import load_delta50_molecules, load_experimental_shifts; ..."
Experimental: 50 molecules
Structures: 50 molecules
All checks passed!
```

## Next Phase Readiness

**Ready for:** Phase 08 Plan 02 (benchmark runner)

**Dependencies satisfied:**
- [x] 50 XYZ files exist and load correctly
- [x] Experimental shifts available with atom assignments
- [x] Python API ready for benchmark execution

**No blockers identified.**
