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

