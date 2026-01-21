# Phase 7: NWChem Integration - Research

**Researched:** 2026-01-21
**Domain:** Computational chemistry - NWChem quantum chemistry package integration
**Confidence:** MEDIUM

## Summary

Replacing ISiCLE wrapper with direct NWChem integration requires understanding four key areas: (1) NWChem input file generation for geometry optimization and NMR shielding calculations, (2) COSMO solvation model syntax, (3) RDKit for SMILES-to-3D conversion and geometry file handling, and (4) parsing NWChem output files for shielding tensors.

**Key findings:**
- NWChem uses a directive-based input format with separate tasks for geometry optimization (`task dft optimize`) and property calculations (`task dft property`)
- COSMO solvation requires a simple `dielec` parameter (4.8 for CHCl3, 46.0 for DMSO) applied to both optimization and NMR steps
- RDKit provides complete SMILES-to-3D pipeline via `EmbedMolecule()` with ETKDG and `UFFOptimizeMolecule()`, plus XYZ/SDF I/O
- NWChem output contains "Chemical Shielding Tensors (GIAO, in ppm)" section with per-atom isotropic values requiring regex parsing
- ISiCLE automation tool (github.com/zaid-shekhani/NMR_Chemical_Shift-Automation) demonstrates complete workflow from SMILES to chemical shifts

**Primary recommendation:** Build direct NWChem I/O using Python string templates for input generation and regex-based output parsing. Use RDKit's ETKDG method for 3D embedding (superior to older methods). Apply COSMO to both geometry optimization AND NMR shielding tasks for consistency.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| RDKit | 2025.09.4+ | SMILES-to-3D, XYZ/SDF I/O | Industry standard cheminformatics toolkit, built-in ETKDG conformer generation |
| NWChem | 7.0.2+ | DFT calculations | Open-source quantum chemistry, COSMO solvation, GIAO NMR shielding |
| Python stdlib | 3.10+ | String templating, regex parsing | No dependencies for file I/O and text processing |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| numpy | Latest | Array operations | If implementing tensor manipulation (not needed for isotropic-only) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| RDKit | Open Babel | Open Babel has XYZ I/O but RDKit has better Python integration and ETKDG method |
| Direct parsing | cclib parser | cclib adds dependency but provides robust multi-package parsing; overkill for NWChem-only |
| String templates | Jinja2 templates | Jinja2 cleaner for complex templates but adds dependency; NWChem inputs are simple enough |

**Installation:**
```bash
# RDKit via conda (recommended)
conda install -c conda-forge rdkit

# Or pip (may have limitations)
pip install rdkit

# NWChem via system package manager
# Ubuntu/Debian: apt install nwchem
# Conda: conda install -c conda-forge nwchem
```

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/
├── nwchem/                 # New module replacing isicle_wrapper.py
│   ├── __init__.py         # Public API
│   ├── input_gen.py        # NWChem input file generation
│   ├── output_parser.py    # Parse NWChem output for shielding
│   ├── geometry.py         # SMILES → 3D, XYZ/SDF handling via RDKit
│   └── templates.py        # Input file string templates
└── [existing modules]
```

**Alternative (single file):**
```
src/qm_nmr_calc/
├── nwchem.py              # All NWChem functionality in one module
└── [existing modules]
```

### Pattern 1: SMILES-to-3D Conversion with RDKit
**What:** Convert SMILES to 3D coordinates using ETKDG method and UFF optimization
**When to use:** All new calculations starting from SMILES
**Example:**
```python
# Source: https://www.rdkit.org/docs/GettingStartedInPython.html
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_3d(smiles: str) -> Chem.Mol:
    """Convert SMILES to 3D molecule with optimized geometry."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Add hydrogens (required for realistic 3D geometry)
    mol = Chem.AddHs(mol)

    # Embed with ETKDG (uses CSD torsion preferences)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xf00d  # Reproducible conformers
    result = AllChem.EmbedMolecule(mol, params)
    if result != 0:
        raise RuntimeError("Failed to embed molecule")

    # UFF force field optimization
    AllChem.UFFOptimizeMolecule(mol)

    return mol
```

### Pattern 2: NWChem Input Generation via String Templates
**What:** Generate NWChem input files using Python f-strings or template strings
**When to use:** All geometry optimization and NMR shielding tasks
**Example:**
```python
# Source: Research synthesis from https://nwchemgit.github.io/Getting-Started.html
def generate_optimization_input(
    geometry_xyz: str,
    functional: str,
    basis_set: str,
    cosmo_dielec: float,
    max_iter: int = 150,
) -> str:
    """Generate NWChem input for geometry optimization with COSMO."""
    return f"""start molecule
title "Geometry Optimization"

geometry units angstrom noautosym
{geometry_xyz}
end

basis spherical
  * library "{basis_set}"
end

dft
  xc {functional}
  iterations {max_iter}
  convergence energy 1e-7 density 1e-5 gradient 5e-4
  direct
end

cosmo
  dielec {cosmo_dielec}
end

task dft optimize
"""
```

### Pattern 3: NWChem Output Parsing with Regex
**What:** Extract isotropic shielding values from NWChem output
**When to use:** After every NMR shielding calculation
**Example:**
```python
# Source: Pattern derived from https://github.com/nwchemgit/nwchem/blob/master/contrib/parsers/nw_spectrum.py
import re
from typing import List, Dict

def parse_shielding_output(output_text: str) -> List[Dict]:
    """
    Parse NWChem output for Chemical Shielding Tensors (GIAO, in ppm).

    Returns list of dicts with keys: atom_index, element, isotropic_ppm
    """
    results = []

    # Find the shielding section
    shielding_section = re.search(
        r'Chemical Shielding Tensors \(GIAO, in ppm\)(.*?)(?=\n\s*\n|\Z)',
        output_text,
        re.DOTALL
    )

    if not shielding_section:
        raise RuntimeError("Shielding section not found in NWChem output")

    # Parse each atom's shielding
    # Pattern matches lines like "Atom:  1  C  Isotropic:  183.4567"
    # Actual NWChem format may vary - needs verification against real output
    atom_pattern = re.compile(
        r'Atom:\s+(\d+)\s+([A-Z][a-z]?)\s+.*?Isotropic:\s+([-\d.]+)',
        re.MULTILINE
    )

    for match in atom_pattern.finditer(shielding_section.group(1)):
        results.append({
            'atom_index': int(match.group(1)),
            'element': match.group(2),
            'isotropic_ppm': float(match.group(3))
        })

    return results
```

### Pattern 4: XYZ/SDF File Handling for Pre-optimized Geometries
**What:** Read pre-optimized geometries from XYZ or SDF files
**When to use:** When `skip_optimization=True` parameter is provided
**Example:**
```python
# Source: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

def load_geometry_from_file(filepath: str, charge: int = 0) -> Chem.Mol:
    """
    Load molecular geometry from XYZ or SDF file.

    Args:
        filepath: Path to .xyz or .sdf file
        charge: Molecular charge (required for XYZ bond determination)
    """
    ext = filepath.lower().split('.')[-1]

    if ext == 'sdf':
        # SDF contains bond info, straightforward
        mol = Chem.MolFromMolFile(filepath, removeHs=False, strictParsing=True)
        if mol is None:
            raise ValueError(f"Failed to parse SDF file: {filepath}")
        return mol

    elif ext == 'xyz':
        # XYZ only has coordinates, must determine bonds
        mol = Chem.MolFromXYZFile(filepath)
        if mol is None:
            raise ValueError(f"Failed to parse XYZ file: {filepath}")

        # Determine bonds from coordinates (requires charge!)
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
        return mol

    else:
        raise ValueError(f"Unsupported file format: {ext}")
```

### Pattern 5: Two-Step DFT Calculation (Optimize + NMR)
**What:** Run geometry optimization, then NMR shielding on optimized geometry
**When to use:** Standard workflow for all SMILES-based calculations
**Example:**
```python
# Source: Research synthesis from current isicle_wrapper.py workflow
def run_two_step_calculation(
    smiles: str,
    job_dir: Path,
    preset: dict,
    solvent_dielec: float,
) -> dict:
    """Run geometry optimization followed by NMR shielding."""

    # Step 1: SMILES to 3D
    mol = smiles_to_3d(smiles)
    xyz_block = Chem.MolToXYZBlock(mol)

    # Step 2: Geometry optimization
    opt_input = generate_optimization_input(
        geometry_xyz=xyz_block,
        functional=preset['functional'],
        basis_set=preset['basis_set'],
        cosmo_dielec=solvent_dielec,
        max_iter=preset['max_iter']
    )

    opt_input_file = job_dir / 'optimize.nw'
    opt_input_file.write_text(opt_input)

    # Run NWChem optimization
    run_nwchem(opt_input_file, job_dir / 'optimize.out')

    # Parse optimized geometry from output
    optimized_xyz = extract_optimized_geometry(job_dir / 'optimize.out')

    # Step 3: NMR shielding on optimized geometry
    nmr_input = generate_shielding_input(
        geometry_xyz=optimized_xyz,
        functional=preset['functional'],
        basis_set=preset['nmr_basis_set'],
        cosmo_dielec=solvent_dielec,
    )

    nmr_input_file = job_dir / 'shielding.nw'
    nmr_input_file.write_text(nmr_input)

    # Run NWChem NMR
    run_nwchem(nmr_input_file, job_dir / 'shielding.out')

    # Parse shielding results
    shielding_data = parse_shielding_output(
        (job_dir / 'shielding.out').read_text()
    )

    return {
        'optimized_geometry': optimized_xyz,
        'shielding_data': shielding_data,
    }
```

### Anti-Patterns to Avoid

- **Mixing gas-phase and COSMO calculations:** Apply COSMO to BOTH optimization and NMR steps, not just one (current bug in isicle_wrapper.py)
- **Forgetting hydrogens in 3D embedding:** Always call `Chem.AddHs()` before `EmbedMolecule()` for realistic geometries
- **Hardcoding basis set names without quotes:** Use `"6-31G*"` not `6-31G*` in NWChem input (special characters require quoting)
- **Inferring skip_optimization from file format:** Require explicit `skip_optimization=True` parameter, don't guess from .xyz/.sdf extension
- **Not validating SMILES before 3D conversion:** Check `Chem.MolFromSmiles()` returns non-None before proceeding
- **Assuming XYZ files include bond orders:** XYZ only has coordinates; use `rdDetermineBonds.DetermineBonds()` or prefer SDF format

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SMILES-to-3D conformers | Custom force field, torsion sampling | RDKit's ETKDG method (`AllChem.ETKDGv3()`) | ETKDG uses Cambridge Structural Database torsion preferences, extensively validated |
| XYZ bond determination | Distance-based heuristics | RDKit's `rdDetermineBonds.DetermineBonds()` | Handles charged molecules, radicals, edge cases from Google Summer of Code 2022 |
| NWChem basis set library | Custom basis set files | NWChem's built-in library via `* library "basis"` | Comprehensive, tested, maintained by NWChem team |
| Molecular file format conversion | Manual string parsing | RDKit's `MolFromXYZFile`, `MolFromMolFile`, `MolToXYZBlock` | Handles edge cases, header variations, whitespace |
| COSMO solvent parameters | Solvent parameter databases | NWChem's built-in COSMO with dielectric constant only | Recent improvements (2026) make simple `dielec` parameter sufficient |

**Key insight:** RDKit and NWChem have solved the hard problems (conformer generation, bond perception, solvation models). Focus on glue code, not reimplementing quantum chemistry or cheminformatics.

## Common Pitfalls

### Pitfall 1: XYZ Files Without Charge Information
**What goes wrong:** `rdDetermineBonds.DetermineBonds()` fails or assigns wrong bonds when molecular charge is unknown
**Why it happens:** XYZ format only stores coordinates, not charge or connectivity. Bond determination algorithm requires formal charge.
**How to avoid:**
- If accepting XYZ files, require charge as separate parameter
- Prefer SDF format (includes bond orders and charge) for pre-optimized geometries
- Validate bond determination succeeded before using molecule
**Warning signs:**
- Error: "Final molecular charge (0) does not match input (-1)"
- Incorrect bond orders (single bonds instead of double bonds in charged species)

### Pitfall 2: COSMO Applied Inconsistently
**What goes wrong:** Geometry optimized in gas phase, NMR calculated in solvent (or vice versa). Results don't match experimental data.
**Why it happens:** Easy to forget COSMO directive in one task but not the other
**How to avoid:**
- Always use same COSMO parameters for both `task dft optimize` and `task dft property`
- Create helper function that adds COSMO block to all generated inputs
- Current codebase has this bug (`cosmo=False` hardcoded in both tasks)
**Warning signs:**
- Chemical shifts ~10 ppm off from experimental (mentioned in isicle_wrapper.py TODO)
- Output shows "Solvent: none" in one calculation but "Solvent: water" in another

### Pitfall 3: NWChem Convergence Failures in Geometry Optimization
**What goes wrong:** Optimization doesn't converge within max_iter, calculation fails
**Why it happens:** Default max iterations (20-40) insufficient for complex molecules, poor initial geometry
**How to avoid:**
- Set `iterations` parameter in DFT block (current presets use 100-150)
- Use better initial guess via RDKit UFF optimization before DFT
- Increase `max_iter` for production calculations vs draft
- Check convergence criteria (energy, density, gradient) are appropriate
**Warning signs:**
- Output shows "geometry did not converge"
- GMAX or GRMS values not decreasing
- Log shows hitting iteration limit

### Pitfall 4: Basis Set Syntax Errors
**What goes wrong:** NWChem fails to parse basis set name, crashes or uses wrong basis
**Why it happens:** Special characters (*,+) in basis set names, capitalization mismatches
**How to avoid:**
- Always quote basis set names: `* library "6-31G*"` not `* library 6-31G*`
- Use lowercase for standard sets: `"6-31g*"` not `"6-31G*"` (NWChem is case-insensitive but lowercase is convention)
- Test basis set syntax with simple molecule before production runs
**Warning signs:**
- NWChem error: "basis set not found"
- Output shows different basis than expected
- Parse error in BASIS directive

### Pitfall 5: Shielding Tensor Parsing Failures
**What goes wrong:** Regex doesn't match output format, returns empty list or partial results
**Why it happens:** NWChem output format varies slightly between versions, GIAO section format assumptions wrong
**How to avoid:**
- Test regex against actual NWChem 7.0.2+ output files first
- Make regex flexible to whitespace variations
- Verify all expected atoms appear in parsed results
- Log warning if atom count doesn't match geometry
**Warning signs:**
- Parsed shielding list shorter than expected
- Missing specific elements (e.g., no hydrogens found)
- Regex matches zero results despite successful calculation

### Pitfall 6: RDKit Embedding Failures for Strained Molecules
**What goes wrong:** `AllChem.EmbedMolecule()` returns -1 (failure), no 3D coordinates generated
**Why it happens:** Molecule too constrained, rings cause conflicts, ETKDG can't find valid geometry
**How to avoid:**
- Check `EmbedMolecule()` return value (0 = success, -1 = failure)
- Try multiple random seeds if first attempt fails
- Consider `useRandomCoords=True` parameter for difficult molecules
- Fall back to simpler embedding without ETKDG as last resort
**Warning signs:**
- `EmbedMolecule()` returns -1
- Large strained rings or fused cage structures
- Molecules with unusual valences

### Pitfall 7: COSMO Solvent Parameter Mistakes
**What goes wrong:** Wrong dielectric constant used, unsupported solvent requested
**Why it happens:** Dielectric constants lookup errors, user requests solvent without known parameters
**How to avoid:**
- Hardcode supported solvents: `SOLVENTS = {'CHCl3': 4.8, 'DMSO': 46.0}`
- Reject unsupported solvents with clear error message listing valid options
- Document dielectric constant sources in comments
- Don't allow custom dielectric constants (prevents user errors)
**Warning signs:**
- User requests "chloroform" but code expects "CHCl3"
- Dielectric constant significantly different from literature (4.8 for CHCl3, 46.0 for DMSO)

## Code Examples

Verified patterns from official sources:

### Complete NWChem Input: Geometry Optimization with COSMO
```nwchem
# Source: https://nwchemgit.github.io/Getting-Started.html
# https://nwchemgit.github.io/Solvation-Models.html
start ethanol_opt
title "Ethanol Geometry Optimization in Chloroform"

geometry units angstrom noautosym
  C        0.00000000     0.00000000     0.00000000
  C        1.52190000     0.00000000     0.00000000
  O        2.04990000     1.31120000     0.00000000
  H       -0.36620000    -0.52390000     0.89000000
  H       -0.36620000    -0.52390000    -0.89000000
  H       -0.36620000     1.02790000     0.00000000
  H        1.88810000    -0.52390000     0.89000000
  H        1.88810000    -0.52390000    -0.89000000
  H        2.72610000     1.31120000     0.00000000
end

basis spherical
  * library "6-31G*"
end

dft
  xc b3lyp
  iterations 150
  convergence energy 1e-7 density 1e-5 gradient 5e-4
  direct
end

cosmo
  dielec 4.8
end

task dft optimize
```

### Complete NWChem Input: NMR Shielding with COSMO
```nwchem
# Source: https://nwchemgit.github.io/Properties.html
# https://nwchemgit.github.io/Solvation-Models.html
start ethanol_nmr
title "Ethanol NMR Shielding in Chloroform"

geometry units angstrom noautosym
  [optimized geometry from previous calculation]
end

basis spherical
  * library "6-311+G(2d,p)"
end

dft
  xc b3lyp
  direct
end

cosmo
  dielec 4.8
end

property
  shielding
end

task dft property
```

### RDKit: Complete SMILES-to-XYZ Workflow
```python
# Source: https://www.rdkit.org/docs/GettingStartedInPython.html
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_xyz(smiles: str, output_file: str) -> None:
    """Convert SMILES to XYZ file with 3D coordinates."""
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates with ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xf00d
    result = AllChem.EmbedMolecule(mol, params)

    if result != 0:
        raise RuntimeError("Failed to embed molecule")

    # Optimize with UFF
    AllChem.UFFOptimizeMolecule(mol)

    # Write XYZ
    xyz_block = Chem.MolToXYZBlock(mol)
    with open(output_file, 'w') as f:
        f.write(xyz_block)
```

### RDKit: Read Pre-optimized Geometry
```python
# Source: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html
# https://rdkit.org/docs/source/rdkit.Chem.rdDetermineBonds.html
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from pathlib import Path

def read_geometry(filepath: Path, charge: int = 0) -> Chem.Mol:
    """
    Read molecular geometry from XYZ or SDF file.

    Args:
        filepath: Path to geometry file
        charge: Molecular formal charge (required for XYZ)

    Returns:
        RDKit Mol object with 3D coordinates
    """
    ext = filepath.suffix.lower()

    if ext == '.sdf':
        mol = Chem.MolFromMolFile(str(filepath), removeHs=False)
        if mol is None:
            raise ValueError(f"Failed to parse SDF: {filepath}")
        return mol

    elif ext == '.xyz':
        mol = Chem.MolFromXYZFile(str(filepath))
        if mol is None:
            raise ValueError(f"Failed to parse XYZ: {filepath}")

        # Determine bonds from coordinates
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=charge)
        except Exception as e:
            raise RuntimeError(f"Bond determination failed: {e}")

        return mol

    else:
        raise ValueError(f"Unsupported format: {ext}")
```

### Python: Extract Optimized Geometry from NWChem Output
```python
# Source: Pattern from https://github.com/nwchemgit/nwchem/blob/master/contrib/parsers/nw_spectrum.py
import re

def extract_optimized_geometry(output_file: str) -> str:
    """
    Extract final optimized geometry from NWChem optimization output.

    Returns XYZ-format string (without element count header).
    """
    with open(output_file) as f:
        content = f.read()

    # Find final geometry section
    # Pattern: "Output coordinates in angstroms" followed by atom lines
    match = re.search(
        r'Output coordinates in angstroms.*?\n\s*\n(.*?)\n\s*\n',
        content,
        re.DOTALL
    )

    if not match:
        raise RuntimeError("Optimized geometry not found in output")

    # Parse atom lines: "  C   0.0000  0.0000  0.0000"
    geom_section = match.group(1)
    xyz_lines = []

    for line in geom_section.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 4:
            element = parts[0]
            x, y, z = parts[1:4]
            xyz_lines.append(f"{element:2s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}")

    return '\n'.join(xyz_lines)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Distance geometry | ETKDG (Experimental Torsion Knowledge Distance Geometry) | 2018.09 (RDKit) | Better conformer quality using CSD torsion preferences |
| ETKDGv1/v2 | ETKDGv3 | 2024.03 (RDKit) | Default method, improved performance |
| Manual XYZ bond parsing | `rdDetermineBonds` module | 2022.09 (RDKit) | Robust bond determination from Google Summer of Code project |
| Basic COSMO | Improved COSMO with new cavity construction | 2026 (NWChem 7.3.0) | Substantially improved solvation energy accuracy |
| Isolated ISiCLE usage | Direct NWChem integration | This phase | Remove runtime dependency, enable COSMO solvation |

**Deprecated/outdated:**
- **ISiCLE as runtime dependency**: Now only a reference during development, fully replaced by direct NWChem I/O
- **RDKit distance geometry without ETKDG**: Older embedding methods inferior to ETKDG for drug-like molecules
- **NWChem COSMO with `iscren` parameter**: Deprecated, use `screen` parameter instead
- **Gas-phase NMR with CHESHIRE scaling**: Phase 8 will derive solvent-specific scaling factors, making COSMO mandatory

## Open Questions

Things that couldn't be fully resolved:

1. **Exact NWChem Output Format for Shielding Tensors**
   - What we know: Section header is "Chemical Shielding Tensors (GIAO, in ppm)", contains isotropic values per atom
   - What's unclear: Exact line format, whitespace patterns, whether atom index/element appears on same line as isotropic value
   - Recommendation: Test regex against actual NWChem 7.0.2 output file in early implementation, adjust patterns as needed. GitHub automation tool (zaid-shekhani/NMR_Chemical_Shift-Automation) has examples in Resources/ folder.

2. **NWChem Optimized Geometry Output Format**
   - What we know: Geometry appears after optimization completes, in angstroms
   - What's unclear: Section header keywords, whether "Output coordinates" or "Final geometry" or other marker
   - Recommendation: Examine actual optimization output file, create robust regex that handles version variations. May need to parse from restart file (.db) instead of text output.

3. **ISiCLE Internal Geometry Object Format**
   - What we know: ISiCLE uses internal geometry representation between optimize and shielding tasks
   - What's unclear: Whether it's RDKit Mol, simple XYZ string, or custom class
   - Recommendation: Not critical - we'll use XYZ strings as intermediate format between our optimization and NMR tasks. Simpler than object passing.

4. **COSMO Performance Impact**
   - What we know: COSMO adds computational cost, recent NWChem improvements (2026) enhanced accuracy
   - What's unclear: Typical runtime increase percentage, memory requirements vs gas-phase
   - Recommendation: Benchmark single test molecule (e.g., ethanol) with/without COSMO. Document typical overhead in user-facing docs.

5. **Charge Parameter Handling for Pre-optimized Geometries**
   - What we know: XYZ bond determination requires molecular charge, SDF doesn't
   - What's unclear: How to get charge from user for XYZ uploads, whether to infer from SMILES if also provided
   - Recommendation: If SMILES provided alongside XYZ, compute charge from SMILES canonical form. Otherwise, require charge as API parameter with default=0.

## Sources

### Primary (HIGH confidence)
- [NWChem Official Documentation - Getting Started](https://nwchemgit.github.io/Getting-Started.html) - Input file structure, basis sets
- [NWChem Official Documentation - Geometry Optimization](https://nwchemgit.github.io/Geometry-Optimization.html) - DFT optimization syntax
- [NWChem Official Documentation - Solvation Models](https://nwchemgit.github.io/Solvation-Models.html) - Complete COSMO syntax and parameters
- [NWChem Official Documentation - Properties](https://nwchemgit.github.io/Properties.html) - NMR shielding calculation syntax
- [NWChem Official Documentation - DFT for Molecules](https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html) - DFT directive details, convergence
- [RDKit Official Documentation - Getting Started in Python](https://www.rdkit.org/docs/GettingStartedInPython.html) - SMILES-to-3D, ETKDG, force fields
- [RDKit Official Documentation - rdmolfiles module](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html) - XYZ/SDF file I/O API
- [RDKit Official Documentation - rdDetermineBonds module](https://rdkit.org/docs/source/rdkit.Chem.rdDetermineBonds.html) - Bond determination from coordinates
- [RDKit Blog - Introducing rdDetermineBonds](https://greglandrum.github.io/rdkit-blog/posts/2022-12-18-introducing-rdDetermineBonds.html) - Implementation details, usage patterns

### Secondary (MEDIUM confidence)
- [Journal of Cheminformatics - ISiCLE NMR Framework](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8) - ISiCLE methodology paper (couldn't access full text due to redirect, abstract only)
- [NWChem GitHub - nw_spectrum.py Parser](https://github.com/nwchemgit/nwchem/blob/master/contrib/parsers/nw_spectrum.py) - Example parsing patterns for NWChem output
- [GitHub - NMR Chemical Shift Automation Tool](https://github.com/zaid-shekhani/NMR_Chemical_Shift-Automation) - SMILES-to-NWChem automation, demonstrates output parsing
- [NWChem FAQ](https://www.cenapad.unicamp.br/parque/manuais/NWChem/NWChem_FAQ.html) - Common errors and solutions
- [Journal of Chemical Theory and Computation - NWChem COSMO Improvements](https://pubs.acs.org/doi/10.1021/acs.jctc.5c01368) - Recent COSMO enhancements (2026)

### Tertiary (LOW confidence)
- WebSearch results for dielectric constants (CHCl3: 4.8, DMSO: 46.0) - Multiple sources agree but not from NWChem official docs
- WebSearch results for NWChem output format patterns - Described but no verbatim examples found
- Community forum discussions about convergence issues - Anecdotal, not systematic documentation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - RDKit and NWChem are official, well-documented projects with stable APIs
- Architecture: MEDIUM - Patterns synthesized from examples, not from ISiCLE source code inspection (couldn't access qm/ module directly)
- Pitfalls: MEDIUM - Based on documentation warnings, FAQ, and GitHub issues; haven't run into these personally
- Code examples: HIGH for RDKit/NWChem official examples, MEDIUM for output parsing (regex needs validation against real output)
- COSMO solvation: HIGH - Official NWChem documentation comprehensive, simple dielec parameter well-documented
- Output parsing: LOW - Could not obtain verbatim NWChem output examples, regex patterns are educated guesses

**Research gaps:**
- Could not access ISiCLE source code directly (GitHub file view returned 404 for qm/ module)
- No complete NWChem shielding output example found (section header known, atom line format unclear)
- ISiCLE paper abstract accessible but full methodology behind paywall/redirect

**Research date:** 2026-01-21
**Valid until:** 30 days (stable domain - NWChem and RDKit change slowly, but verify output format assumptions early)
