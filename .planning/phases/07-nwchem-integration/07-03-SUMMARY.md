---
phase: 07-nwchem-integration
plan: 03
subsystem: nwchem/geometry
tags: [rdkit, smiles, xyz, sdf, geometry, conformer]

dependency_graph:
  requires: []
  provides:
    - smiles_to_xyz: SMILES to 3D conversion with ETKDGv3
    - load_geometry_file: XYZ/SDF file loading with bond determination
    - mol_to_xyz_block: NWChem-ready coordinate block output
    - validate_geometry: 3D conformer validation
  affects:
    - 07-04: Runner will use geometry functions for molecule preparation
    - 07-05: Integration will use geometry module for end-to-end workflow

tech_stack:
  added: []
  patterns:
    - ETKDGv3 conformer generation with deterministic seeding
    - UFF force field optimization for initial geometry
    - rdDetermineBonds for XYZ bond determination with charge parameter

key_files:
  created:
    - src/qm_nmr_calc/nwchem/geometry.py
    - tests/test_nwchem_geometry.py
    - tests/fixtures/ethanol.xyz
    - tests/fixtures/ethanol.sdf
  modified:
    - src/qm_nmr_calc/nwchem/__init__.py

decisions:
  - id: etkdgv3-default
    choice: Use ETKDGv3 as default conformer generation method
    reason: Best available method per RDKit docs, uses CSD torsion preferences
  - id: deterministic-seed
    choice: Default random_seed=0xF00D for reproducible conformers
    reason: Enables reproducible calculations while allowing override
  - id: xyz-charge-param
    choice: XYZ bond determination requires explicit charge parameter (default=0)
    reason: rdDetermineBonds algorithm needs formal charge for correct bond orders

metrics:
  duration: 3 min
  completed: 2026-01-21
---

# Phase 07 Plan 03: Geometry Handling Module Summary

RDKit-based geometry preparation for NWChem calculations supporting SMILES-to-3D conversion and pre-optimized geometry loading.

## What Was Built

### Geometry Handling Module (`src/qm_nmr_calc/nwchem/geometry.py`)

**smiles_to_xyz(smiles: str, random_seed: int = 0xF00D) -> tuple[Mol, str]**
- Parses SMILES string and validates
- Adds hydrogens with `Chem.AddHs()` for realistic 3D geometry
- Generates 3D coordinates using ETKDGv3 (best available method)
- Runs UFF force field optimization
- Returns RDKit Mol object and XYZ coordinate block for NWChem

**load_geometry_file(filepath: Path | str, charge: int = 0) -> tuple[Mol, str]**
- Detects format from file extension (.xyz or .sdf)
- For SDF: Loads with bond information preserved
- For XYZ: Loads and determines bonds using rdDetermineBonds with charge parameter
- Validates geometry before returning
- Returns RDKit Mol object and XYZ coordinate block

**mol_to_xyz_block(mol: Mol) -> str**
- Converts RDKit mol to NWChem geometry block format
- Strips first two lines (atom count + title) from standard XYZ
- Returns just "Element X Y Z" lines suitable for NWChem input

**validate_geometry(mol: Mol) -> bool**
- Ensures molecule has at least one conformer
- Verifies coordinates are not all zeros (failed embedding indicator)
- Raises ValueError with clear message if invalid

### Test Coverage

18 unit tests covering:
- SMILES to XYZ conversion (ethanol, methane)
- Hydrogen addition verification
- Reproducibility with same random seed
- Invalid SMILES error handling
- XYZ file loading with bond determination
- SDF file loading with bond preservation
- Unsupported format rejection
- Output format compatibility with NWChem

### Test Fixtures

- `tests/fixtures/ethanol.xyz`: 9-atom XYZ file for testing
- `tests/fixtures/ethanol.sdf`: Same molecule with bond information for SDF loading tests

## Technical Details

### ETKDGv3 Selection

ETKDGv3 (Experimental Torsion Knowledge Distance Geometry v3) was chosen because:
- Uses Cambridge Structural Database torsion angle preferences
- Produces better initial conformers than older distance geometry methods
- Default in RDKit 2024.03+, well-validated

### XYZ Bond Determination

For XYZ files (which contain only coordinates, no connectivity):
- Uses `rdDetermineBonds.DetermineBonds()` from RDKit
- Requires molecular formal charge parameter for correct bond orders
- Default charge=0 for neutral molecules
- Users must specify charge for ions/charged species

### Output Format

XYZ coordinate block format for NWChem:
```
C     0.000000    0.000000    0.000000
C     1.521900    0.000000    0.000000
O     2.049900    1.311200    0.000000
...
```

No atom count or title line - just element + 3 Cartesian coordinates per line.

## Commits

| Hash | Type | Description |
|------|------|-------------|
| 0d610be | feat | Create geometry handling module |
| e61d7bf | test | Add unit tests for geometry handling |
| a8c117e | feat | Export geometry functions from nwchem package |

## Deviations from Plan

None - plan executed exactly as written.

## Next Phase Readiness

**Dependencies satisfied:**
- smiles_to_xyz provides SMILES-to-3D conversion for new submissions
- load_geometry_file supports pre-optimized XYZ/SDF uploads
- Output format matches NWChem geometry directive requirements

**Ready for:**
- Plan 04 (Runner) can use geometry functions for molecule preparation
- Plan 05 (Integration) can assemble full workflow with geometry module
