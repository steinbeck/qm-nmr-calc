# Library Documentation

This guide documents how third-party Python libraries are integrated into QM NMR Calculator. Code examples are extracted from the actual codebase, not generic library usage patterns.

**Target audience:** Developers and contributors working on QM NMR Calculator internals.

**Documentation approach:**
- Each section shows WHERE the library is used (integration points table)
- Code examples come from actual source files with line references
- Links to official library docs for comprehensive API reference

**Libraries covered:**
- [RDKit](#rdkit) - SMILES parsing, conformer generation, visualization
- [NWChem](#nwchem) - DFT calculations and NMR shielding
- [Huey](#huey-task-queue) - Async job queue with SQLite persistence

For JavaScript libraries (3Dmol.js, SmilesDrawer) and optional tools (CREST/xTB), see Phase 29-02.

---

## RDKit

RDKit handles all cheminformatics operations: SMILES parsing, 3D conformer generation, and 2D structure visualization.

**Official docs:** [RDKit Documentation](https://www.rdkit.org/docs/)

### Integration Points

| Module | Functions | Purpose |
|--------|-----------|---------|
| `validation.py` | `validate_smiles()`, `validate_mol_file()` | Input validation |
| `conformers/generator.py` | `generate_conformers_kdg()`, `optimize_conformers_mmff()` | 3D structure generation |
| `visualization.py` | `generate_annotated_structure()` | 2D depiction with shift labels |
| `nwchem/geometry.py` | `smiles_to_xyz()`, `mol_to_xyz_block()` | XYZ coordinate conversion |

### SMILES Validation

The `validate_smiles()` function parses SMILES strings and returns either a valid molecule or an error message.

```python
# Source: src/qm_nmr_calc/validation.py, lines 10-39
from rdkit import Chem, RDLogger

def validate_smiles(smiles: str) -> tuple[Optional[Chem.Mol], Optional[str]]:
    """Validate SMILES and return molecule or error message."""
    # Capture RDKit error messages
    stderr_capture = io.StringIO()
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Temporarily redirect stderr to capture RDKit messages
    old_stderr = sys.stderr
    sys.stderr = stderr_capture

    try:
        mol = Chem.MolFromSmiles(smiles)
    finally:
        sys.stderr = old_stderr

    if mol is None:
        error_msg = stderr_capture.getvalue().strip()
        if not error_msg:
            error_msg = "Invalid SMILES string"
        return None, error_msg

    return mol, None
```

**Key pattern:** RDKit logs errors to stderr rather than raising exceptions. The code captures stderr to provide meaningful error messages to users instead of generic "Invalid SMILES" feedback.

### Conformer Generation (KDG)

Conformer generation uses RDKit's Knowledge-based Distance Geometry (KDG) algorithm.

**Why KDG instead of ETKDG?**
- **ETKDG** (Experimental-Torsion-angle preference with Distance Geometry) includes torsion preferences derived from crystal structures
- **KDG** uses only distance geometry without crystal structure bias
- For solution-phase NMR predictions, avoiding crystal packing artifacts is preferred

```python
# Source: src/qm_nmr_calc/conformers/generator.py, lines 51-105
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

def generate_conformers_kdg(
    smiles: str, num_confs: int, random_seed: int = 0xF00D
) -> Chem.Mol:
    """Generate conformers using KDG (no crystal bias)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    # Configure KDG parameters
    params = rdDistGeom.KDG()
    params.randomSeed = random_seed
    params.numThreads = 0  # Use all threads
    params.pruneRmsThresh = -1.0  # No pre-optimization pruning

    conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    # Fallback: retry with random coordinates for difficult molecules
    if len(conf_ids) == 0:
        params.useRandomCoords = True
        conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

        if len(conf_ids) == 0:
            raise RuntimeError(
                f"Failed to generate conformers for {smiles} even with random coordinates"
            )

    return mol
```

**Key parameters:**
- `randomSeed=0xF00D` - Project-standard seed for reproducibility
- `numThreads=0` - Use all available CPU threads
- `pruneRmsThresh=-1.0` - Disable pre-optimization RMSD pruning (pruning happens after MMFF optimization)
- `useRandomCoords=True` - Fallback for molecules that fail standard embedding

### 2D Structure Drawing

The `generate_annotated_structure()` function creates structure images with chemical shift labels on each atom.

```python
# Source: src/qm_nmr_calc/visualization.py, lines 71-119
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def generate_annotated_structure(
    smiles: str,
    h1_shifts: list[dict],
    c13_shifts: list[dict],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Generate structure with chemical shift annotations."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # CRITICAL: Set atomNote BEFORE calling PrepareMolForDrawing
    # Convert NWChem 1-based indices to RDKit 0-based
    for shift_data in h1_shifts + c13_shifts:
        rdkit_idx = shift_data["index"] - 1
        atom = mol.GetAtomWithIdx(rdkit_idx)
        atom.SetProp("atomNote", f"{shift_data['shift']:.2f}")

    # Prepare molecule for drawing (after setting notes)
    rdMolDraw2D.PrepareMolForDrawing(mol)

    # Generate SVG
    d_svg = rdMolDraw2D.MolDraw2DSVG(800, 600)
    d_svg.drawOptions().annotationFontScale = 0.9
    d_svg.DrawMolecule(mol)
    d_svg.FinishDrawing()

    svg_path = output_dir / "structure_annotated.svg"
    svg_path.write_text(d_svg.GetDrawingText())

    return svg_path, png_path
```

**Critical ordering:** The `atomNote` property MUST be set BEFORE calling `PrepareMolForDrawing()`. Setting it after will have no effect on the rendered output.

**Index conversion:** NWChem uses 1-based atom indices; RDKit uses 0-based. Always subtract 1 when converting.

---

## NWChem

NWChem performs the quantum chemistry calculations: DFT geometry optimization and GIAO NMR shielding tensor computation.

**Official docs:** [NWChem Documentation](https://nwchemgit.github.io/)

### Integration Points

| Module | Functions | Purpose |
|--------|-----------|---------|
| `nwchem/input_gen.py` | `generate_optimization_input()`, `generate_shielding_input()` | Input file creation |
| `nwchem/output_parser.py` | `parse_shielding_output()`, `extract_dft_energy()`, `extract_optimized_geometry()` | Result extraction |
| `nwchem/runner.py` | `run_nwchem()` | Subprocess execution |
| `nwchem/geometry.py` | `smiles_to_xyz()`, `mol_to_xyz_block()` | XYZ coordinate conversion |

### Input File Generation

Input files are generated using a template-based approach with f-strings. The format follows ISiCLE conventions for reliability.

```python
# Source: src/qm_nmr_calc/nwchem/input_gen.py, lines 33-97
def generate_optimization_input(
    geometry_xyz: str,
    functional: str,
    basis_set: str,
    solvent: str,
    max_iter: int = 150,
) -> str:
    """Generate NWChem input file for geometry optimization."""
    solvent_name = _validate_solvent(solvent)

    # Build COSMO block only for non-vacuum calculations
    if solvent_name == "vacuum":
        cosmo_block = ""
    else:
        cosmo_block = f"""
cosmo
  do_gasphase False
  solvent {solvent_name}
end
"""

    return f"""start molecule
title "Geometry Optimization"

memory global 1600 mb heap 100 mb stack 600 mb

geometry units angstrom noautoz noautosym
{geometry_xyz}
end

basis spherical
  * library {basis_set}
end

dft
  xc {functional}
end

driver
  maxiter {max_iter}
  xyz molecule_geom
end
{cosmo_block}
task dft optimize
"""
```

**Key input sections:**
- `memory` - Memory allocation (1600 MB global, 100 MB heap, 600 MB stack)
- `geometry` - XYZ coordinates with `noautoz noautosym` to prevent automatic symmetry detection
- `basis` - Basis set library (e.g., `6-31G*`, `6-311+G(2d,p)`)
- `dft` - DFT functional specification (e.g., `b3lyp`)
- `driver` - Optimization control with iteration limit
- `cosmo` - Implicit solvation block (omitted for vacuum/gas-phase)

**COSMO solvation:** The COSMO block is only included for non-vacuum solvents. Supported solvents: `chcl3`, `dmso`, `water`, `acetone`, `methanol`, `vacuum`.

### Output Parsing

NWChem output is parsed using regex patterns to extract shielding tensors and optimized geometries.

```python
# Source: src/qm_nmr_calc/nwchem/output_parser.py, lines 150-240
import re

def parse_shielding_output(output_text: str) -> dict:
    """Parse NMR shielding tensors from NWChem output.

    Returns:
        {
            "index": [1, 2, 3, ...],      # 1-based atom indices
            "atom": ["C", "H", "H", ...], # Element symbols
            "shielding": [183.4, 29.1, ...]  # Isotropic shielding in ppm
        }
    """
    # Pattern: "Atom:    1  C" followed by "isotropic  =     183.4567"
    atom_block_pattern = re.compile(
        r"Atom:\s+(\d+)\s+([A-Za-z][a-z]?)"  # "Atom:    1  C"
        r".*?"                                # Any content between
        r"isotropic\s*=\s*([-\d.]+)",         # "isotropic  =     183.4567"
        re.DOTALL | re.IGNORECASE,
    )

    indices, atoms, shieldings = [], [], []
    for match in atom_block_pattern.finditer(output_text):
        indices.append(int(match.group(1)))
        atoms.append(match.group(2))
        shieldings.append(float(match.group(3)))

    return {"index": indices, "atom": atoms, "shielding": shieldings}
```

**NWChem output format:**
```
GIAO Chemical Shielding Tensors
--------------------------------
Atom:    1  C
   isotropic  =     183.4567
   ...
```

**Key patterns:**
- Atom block starts with `Atom:` followed by 1-based index and element symbol
- Isotropic shielding appears on a subsequent line as `isotropic = value`
- The `re.DOTALL` flag allows `.` to match newlines between atom header and isotropic value

### Why Direct NWChem Integration

The project originally used ISiCLE (a Python wrapper for NWChem) in v1.0 but switched to direct NWChem calls in v2.0:

**Reasons for direct integration:**
1. **Better error handling** - Direct access to NWChem stdout/stderr for debugging
2. **Fine-grained control** - Full control over input file parameters and memory allocation
3. **Clearer debugging** - Raw `.out` files preserved in scratch directories for inspection
4. **Simpler dependencies** - No need to maintain ISiCLE compatibility

**Scratch directory structure:** Each job writes NWChem files to `./data/scratch/{job_id}/`:
- `optimize.nw` - Geometry optimization input
- `optimize.out` - Optimization output (preserved for debugging)
- `shielding.nw` - NMR shielding input
- `shielding.out` - Shielding output with tensor data

See [Architecture Documentation](architecture.md#scratch-directory) for detailed scratch directory layout.

