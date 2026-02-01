# Phase 29: Library Documentation - Research

**Researched:** 2026-02-01
**Domain:** Technical Documentation / Library Integration Guides
**Confidence:** HIGH

## Summary

This phase documents how six third-party libraries are integrated into the QM NMR Calculator codebase. The research involved deep analysis of the existing source code to understand actual usage patterns, combined with documentation best practices research.

The codebase already has mature, well-structured integrations for all six libraries:
- **RDKit** (~4 modules): SMILES parsing, conformer generation, visualization
- **NWChem** (~4 modules): Input generation, subprocess execution, output parsing
- **Huey** (~2 modules): Task queue with signal-based status updates
- **3Dmol.js** (~2 templates): Interactive 3D viewer with shift labels
- **CREST/xTB** (~1 module): Optional ensemble conformer generation
- **SmilesDrawer** (~1 template): Real-time SMILES preview

The documentation task is primarily a **documentation extraction exercise** - the integration patterns already exist and work correctly. The goal is to make these patterns visible to contributors.

**Primary recommendation:** Use the Diataxis framework (tutorials/how-to/reference/explanation) to structure each library section, with code examples extracted directly from the existing codebase.

## Standard Stack

The established documentation approach for this domain:

### Core
| Tool | Purpose | Why Standard |
|------|---------|--------------|
| Markdown | Documentation format | Universal, GitHub-rendered, cross-platform |
| Diataxis Framework | Content organization | Industry-standard 4-quadrant approach |
| Code fences | Example formatting | Syntax highlighting, copy-paste ready |
| Mermaid diagrams | Visual flows | Already used in architecture.md |

### Supporting
| Tool | Purpose | When to Use |
|------|---------|-------------|
| Inline links | Cross-referencing | Link to architecture.md, usage.md for context |
| Tables | API summaries | Function signatures, parameter docs |
| Details/summary | Optional content | Extended examples, edge cases |

### Documentation Structure
Based on the existing docs/ directory:

```
docs/
  libraries.md          # Main library documentation file
```

The docs/README.md already links to libraries.md as "Integration guides for RDKit, NWChem, Huey, and other dependencies" for the developer audience.

## Architecture Patterns

### Recommended Document Structure

Based on the Diataxis framework and existing codebase:

```markdown
# Library Documentation

## RDKit
### Overview (Explanation)
### Core Functions (Reference)
### Usage Examples (How-To)

## NWChem
### Overview (Explanation)
### Input Generation (Reference)
### Output Parsing (Reference)
### Subprocess Execution (Reference)

## Huey Task Queue
### Overview (Explanation)
### Configuration (Reference)
### Signal Handlers (Reference)

## 3Dmol.js
### Overview (Explanation)
### Viewer Setup (Reference)
### Shift Labels (How-To)

## CREST/xTB (Optional)
### Overview (Explanation)
### Availability Detection (Reference)
### Ensemble Parsing (Reference)

## SmilesDrawer
### Overview (Explanation)
### Canvas Setup (Reference)
```

### Pattern 1: Code-First Documentation

**What:** Extract actual code patterns from the codebase, not theoretical examples
**When to use:** All library documentation sections
**Why:** Contributors can trust examples work because they ARE the working code

Example structure:
```markdown
### SMILES Validation

The `validate_smiles()` function in `validation.py` handles SMILES parsing:

```python
# From src/qm_nmr_calc/validation.py
from rdkit import Chem, RDLogger

def validate_smiles(smiles: str) -> tuple[Optional[Chem.Mol], Optional[str]]:
    """Validate SMILES and return molecule or error message."""
    # Capture RDKit error messages
    RDLogger.logger().setLevel(RDLogger.ERROR)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    return mol, None
```
```

### Pattern 2: Integration Point Documentation

**What:** Document WHERE libraries are used, not just HOW
**When to use:** Each library overview section
**Example:**

```markdown
## RDKit

RDKit is used in the following modules:

| Module | Functions | Purpose |
|--------|-----------|---------|
| `validation.py` | `validate_smiles()`, `validate_mol_file()` | Input validation |
| `conformers/generator.py` | `generate_conformers_kdg()`, `optimize_conformers_mmff()` | 3D structure |
| `visualization.py` | `generate_annotated_structure()` | 2D depiction |
| `nwchem/geometry.py` | `smiles_to_xyz()`, `mol_to_xyz_block()` | XYZ conversion |
```

### Pattern 3: External vs Internal APIs

**What:** Clearly distinguish library APIs from our wrapper functions
**When to use:** Reference sections
**Example:**

```markdown
### Conformer Generation

**RDKit API:**
- `rdDistGeom.KDG()` - Knowledge-based Distance Geometry parameters
- `rdDistGeom.EmbedMultipleConfs()` - Multi-conformer embedding

**Our Wrapper:**
- `generate_conformers_kdg(smiles, num_confs, random_seed)` - Project-standard interface
```

### Anti-Patterns to Avoid
- **Generic documentation:** "RDKit can parse SMILES" - instead show exact usage
- **Outdated examples:** Don't invent examples, use actual codebase code
- **Missing context:** Always link to the source file location
- **API duplication:** Don't re-document library APIs, link to official docs

## Don't Hand-Roll

This phase is documentation, not implementation. The following already exist:

| Need | Don't Build | Already Have |
|------|-------------|--------------|
| RDKit patterns | New examples | `validation.py`, `generator.py`, `visualization.py` |
| NWChem patterns | New input formats | `nwchem/input_gen.py`, `output_parser.py` |
| Huey patterns | New task code | `queue.py`, `tasks.py` |
| 3Dmol.js patterns | New viewer code | `results.html`, `status.html` |
| CREST patterns | New subprocess code | `crest_generator.py` |
| SmilesDrawer patterns | New preview code | `submit.html` |

**Key insight:** The documentation task extracts and explains existing patterns, not creates new ones.

## Common Pitfalls

### Pitfall 1: Documenting Library APIs Instead of Integration Patterns

**What goes wrong:** Documentation describes what RDKit/NWChem CAN do instead of what THIS PROJECT does with them
**Why it happens:** Temptation to be comprehensive about the library
**How to avoid:** Focus on "how we use it" not "what it offers"
**Warning signs:** Sections that could apply to any project using the library

### Pitfall 2: Code Examples Drift from Actual Code

**What goes wrong:** Documentation examples diverge from actual implementation over time
**Why it happens:** Examples written separately, not extracted from source
**How to avoid:** Include file paths in every code block, note line numbers if helpful
**Warning signs:** Examples that don't quite match codebase patterns

### Pitfall 3: Missing "Why" Behind Decisions

**What goes wrong:** Documentation says "we use X" but not why X over alternatives
**Why it happens:** Focus on "what" without "why"
**How to avoid:** Include brief rationale for non-obvious choices
**Warning signs:** Contributors asking "why not Y?" for documented patterns

Example of good "why" documentation:
```markdown
### Why KDG Instead of ETKDG?

We use `rdDistGeom.KDG()` instead of `ETKDGv3()` for conformer generation:

- **ETKDG** includes crystal structure bias (torsion preferences from X-ray data)
- **KDG** uses only knowledge-based distance geometry without this bias
- For solution-phase NMR predictions, avoiding crystal packing artifacts is preferred
```

### Pitfall 4: JavaScript Documentation Without Context

**What goes wrong:** 3Dmol.js/SmilesDrawer examples isolated from template context
**Why it happens:** Treating JS libraries same as Python modules
**How to avoid:** Show complete script blocks including CDN includes and DOM setup
**Warning signs:** Code that wouldn't run if copy-pasted into a new template

## Code Examples

All examples below are from the actual codebase. These patterns should be documented.

### RDKit: SMILES Validation (validation.py)

```python
# Source: src/qm_nmr_calc/validation.py, lines 10-39
from rdkit import Chem, RDLogger

def validate_smiles(smiles: str) -> tuple[Optional[Chem.Mol], Optional[str]]:
    """Validate SMILES and return molecule or error message."""
    # Suppress RDKit warnings, capture errors
    RDLogger.logger().setLevel(RDLogger.ERROR)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    return mol, None
```

### RDKit: Conformer Generation (conformers/generator.py)

```python
# Source: src/qm_nmr_calc/conformers/generator.py, lines 51-106
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

def generate_conformers_kdg(smiles: str, num_confs: int, random_seed: int = 0xF00D) -> Chem.Mol:
    """Generate conformers using KDG (no crystal bias)."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # KDG parameters (not ETKDG - no crystal structure bias)
    params = rdDistGeom.KDG()
    params.randomSeed = random_seed
    params.numThreads = 0  # All threads
    params.pruneRmsThresh = -1.0  # No pre-pruning

    conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    if len(conf_ids) == 0:
        params.useRandomCoords = True  # Fallback
        conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    return mol
```

### RDKit: 2D Structure Drawing (visualization.py)

```python
# Source: src/qm_nmr_calc/visualization.py, lines 71-119
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def generate_annotated_structure(smiles: str, h1_shifts: list, c13_shifts: list, output_dir: Path):
    """Generate structure with chemical shift annotations."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Set atomNote BEFORE PrepareMolForDrawing
    for shift_data in h1_shifts + c13_shifts:
        rdkit_idx = shift_data["index"] - 1  # NWChem 1-based -> RDKit 0-based
        atom = mol.GetAtomWithIdx(rdkit_idx)
        atom.SetProp("atomNote", f"{shift_data['shift']:.2f}")

    rdMolDraw2D.PrepareMolForDrawing(mol)

    # SVG output
    d_svg = rdMolDraw2D.MolDraw2DSVG(800, 600)
    d_svg.drawOptions().annotationFontScale = 0.9
    d_svg.DrawMolecule(mol)
    d_svg.FinishDrawing()
```

### NWChem: Input File Generation (nwchem/input_gen.py)

```python
# Source: src/qm_nmr_calc/nwchem/input_gen.py, lines 33-97
def generate_optimization_input(geometry_xyz: str, functional: str, basis_set: str, solvent: str, max_iter: int = 150) -> str:
    """Generate NWChem input for geometry optimization with COSMO."""
    cosmo_block = "" if solvent == "vacuum" else f"""
cosmo
  do_gasphase False
  solvent {solvent}
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

### NWChem: Output Parsing (nwchem/output_parser.py)

```python
# Source: src/qm_nmr_calc/nwchem/output_parser.py, lines 150-240
import re

def parse_shielding_output(output_text: str) -> dict:
    """Parse NMR shielding tensors from NWChem output."""
    # Pattern: "Atom:    1  C" followed by "isotropic  =     183.4567"
    atom_block_pattern = re.compile(
        r"Atom:\s+(\d+)\s+([A-Za-z][a-z]?)"
        r".*?"
        r"isotropic\s*=\s*([-\d.]+)",
        re.DOTALL | re.IGNORECASE,
    )

    indices, atoms, shieldings = [], [], []
    for match in atom_block_pattern.finditer(output_text):
        indices.append(int(match.group(1)))
        atoms.append(match.group(2))
        shieldings.append(float(match.group(3)))

    return {"index": indices, "atom": atoms, "shielding": shieldings}
```

### Huey: Task Queue Configuration (queue.py)

```python
# Source: src/qm_nmr_calc/queue.py, lines 1-106
from huey import SqliteHuey, signals

# SQLite-backed queue with fsync for durability
huey = SqliteHuey('qm-nmr-calc', filename='./data/huey.db', fsync=True)

@huey.signal(signals.SIGNAL_EXECUTING)
def on_task_start(signal, task):
    """Update job status when task starts."""
    job_id = task.args[0]
    update_job_status(job_id, status='running', started_at=datetime.utcnow())

@huey.signal(signals.SIGNAL_COMPLETE)
def on_task_complete(signal, task):
    """Update job status and send notification on success."""
    job_id = task.args[0]
    update_job_status(job_id, status='complete', completed_at=datetime.utcnow())
    # Send email notification if configured

@huey.signal(signals.SIGNAL_ERROR)
def on_task_error(signal, task, exc=None):
    """Update job status with error details on failure."""
    job_id = task.args[0]
    update_job_status(job_id, status='failed', error_message=str(exc))
```

### Huey: Task Definition (tasks.py)

```python
# Source: src/qm_nmr_calc/tasks.py, lines 82-144
from .queue import huey

@huey.task()
def run_nmr_task(job_id: str) -> dict:
    """Execute NMR calculation for a queued job."""
    job_status = load_job_status(job_id)
    preset = PRESETS[PresetName(job_status.input.preset)]

    # Step tracking via callbacks
    start_step(job_id, "geometry_optimization", "Optimizing geometry")

    result = run_calculation(
        smiles=job_status.input.smiles,
        job_dir=get_job_dir(job_id),
        preset=preset,
        solvent=job_status.input.solvent,
        on_optimization_complete=lambda: start_step(job_id, "nmr_shielding", "Computing NMR")
    )

    # Convert shielding to shifts, generate visualizations
    # ...
    return {'success': True, 'job_id': job_id}
```

### 3Dmol.js: Viewer Setup (results.html)

```html
<!-- Source: src/qm_nmr_calc/api/templates/results.html, lines 6-7, 401-507 -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.3/3Dmol-min.js"></script>

<script>
// Create viewer in container
viewer = $3Dmol.createViewer(document.getElementById('viewer-container-3d'), {
    backgroundColor: 'white'
});

// Load geometry (prefer SDF for proper bonds)
const data = await fetch('/api/v1/jobs/' + JOB_ID + '/geometry.json').then(r => r.json());
const model = data.sdf
    ? viewer.addModel(data.sdf, 'sdf')
    : viewer.addModel(data.xyz, 'xyz', {assignBonds: true});

// Apply stick+sphere style
viewer.setStyle({}, {
    stick: {radius: 0.12, colorscheme: 'Jmol'},
    sphere: {scale: 0.25, colorscheme: 'Jmol'}
});

viewer.zoomTo();
viewer.render();
</script>
```

### 3Dmol.js: Shift Labels (results.html)

```javascript
// Source: src/qm_nmr_calc/api/templates/results.html, lines 509-538
function addShiftLabels(model, h1Assignments, c13Assignments) {
    const atoms = model.selectedAtoms({});

    atoms.forEach(atom => {
        // CRITICAL: 3Dmol.js uses 0-based serial, NWChem uses 1-based indices
        const atomIndex = String(atom.serial + 1);

        if (atom.elem === 'H' && h1Assignments && h1Assignments[atomIndex]) {
            viewer.addLabel(h1Assignments[atomIndex].toFixed(2), {
                position: {x: atom.x, y: atom.y, z: atom.z},
                fontSize: 11,
                fontColor: 'white',
                backgroundColor: '#3498db',  // Blue for 1H
                backgroundOpacity: 0.85,
                inFront: true
            });
        }

        if (atom.elem === 'C' && c13Assignments && c13Assignments[atomIndex]) {
            viewer.addLabel(c13Assignments[atomIndex].toFixed(2), {
                position: {x: atom.x, y: atom.y, z: atom.z},
                fontSize: 11,
                fontColor: 'white',
                backgroundColor: '#e67e22',  // Orange for 13C
                backgroundOpacity: 0.85,
                inFront: true
            });
        }
    });
}
```

### CREST/xTB: Availability Detection (crest_generator.py)

```python
# Source: src/qm_nmr_calc/conformers/crest_generator.py, lines 47-64
import shutil
from functools import lru_cache

@lru_cache(maxsize=1)
def detect_crest_available() -> bool:
    """Detect if CREST and xTB binaries are available on PATH."""
    crest_found = shutil.which("crest") is not None
    xtb_found = shutil.which("xtb") is not None
    return crest_found and xtb_found
```

### CREST/xTB: Ensemble Parsing (crest_generator.py)

```python
# Source: src/qm_nmr_calc/conformers/crest_generator.py, lines 98-169
def parse_crest_ensemble(ensemble_file: Path) -> list[CRESTConformer]:
    """Parse CREST multi-structure XYZ ensemble file.

    CREST outputs concatenated XYZ: atom count, energy comment, coordinates.
    """
    conformers = []
    lines = ensemble_file.read_text().strip().split("\n")

    i = 0
    conf_number = 1
    while i < len(lines):
        num_atoms = int(lines[i].strip())
        i += 1
        energy_hartree = float(lines[i].strip().split()[0])  # First token is energy
        i += 1
        coord_lines = [lines[i+j] for j in range(num_atoms)]
        i += num_atoms

        conformer = CRESTConformer(
            conformer_id=f"conf_{conf_number:03d}",
            energy_hartree=energy_hartree,
            xyz_block=f"{num_atoms}\n{energy_hartree}\n" + "\n".join(coord_lines)
        )
        conformers.append(conformer)
        conf_number += 1

    return conformers
```

### SmilesDrawer: Canvas Setup (submit.html)

```html
<!-- Source: src/qm_nmr_calc/api/templates/submit.html, lines 147-214 -->
<script src="https://unpkg.com/smiles-drawer@1.2.0/dist/smiles-drawer.min.js"></script>

<script>
// Configure drawer with custom colors
const smilesDrawer = new SmilesDrawer.Drawer({
    width: 350,
    height: 350,
    bondThickness: 0.6,
    bondLength: 15,
    terminalCarbons: true,
    explicitHydrogens: false,
    themes: {
        light: {
            C: '#222', O: '#e74c3c', N: '#3498db',
            S: '#f1c40f', P: '#e67e22', F: '#1abc9c',
            CL: '#1abc9c', BR: '#e74c3c', I: '#9b59b6', H: '#aaa'
        }
    }
});

// Parse and draw on input change (debounced)
function updateMoleculePreview(smiles) {
    SmilesDrawer.parse(smiles, function(tree) {
        smilesDrawer.draw(tree, 'molecule-preview', 'light', false);
        status.textContent = 'Valid SMILES';
    }, function(err) {
        status.textContent = 'Invalid SMILES syntax';
    });
}
</script>
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| ETKDGv3 for conformers | KDG for NMR | Project design | Avoids crystal packing bias |
| ISiCLE wrapper | Direct NWChem calls | v2.0 | Better error handling |
| Celery + Redis | Huey + SQLite | Project design | Simpler deployment |

**Current library versions (from pyproject.toml):**
- RDKit: >=2025.9.3
- Huey: >=2.6.0
- FastAPI: >=0.128.0

**JavaScript CDN versions (from templates):**
- 3Dmol.js: 2.5.3
- SmilesDrawer: 1.2.0

## Open Questions

None significant. The integration patterns are mature and well-tested. The documentation task is extraction, not research.

Minor consideration:
1. **Level of NWChem detail:** Should the docs explain NWChem input format syntax (e.g., what `noautoz noautosym` means), or just show our templates? Recommendation: Light explanation with links to NWChem official docs.

## Sources

### Primary (HIGH confidence)
- Codebase analysis: Direct reading of all source files
  - `src/qm_nmr_calc/validation.py` - RDKit SMILES/MOL validation
  - `src/qm_nmr_calc/conformers/generator.py` - RDKit conformer generation
  - `src/qm_nmr_calc/visualization.py` - RDKit 2D drawing
  - `src/qm_nmr_calc/nwchem/*.py` - NWChem integration (4 modules)
  - `src/qm_nmr_calc/queue.py` - Huey configuration
  - `src/qm_nmr_calc/tasks.py` - Huey task definitions
  - `src/qm_nmr_calc/conformers/crest_generator.py` - CREST/xTB integration
  - `src/qm_nmr_calc/api/templates/*.html` - 3Dmol.js/SmilesDrawer usage

### Secondary (MEDIUM confidence)
- [Diataxis Framework](https://diataxis.fr/) - Documentation structure
- [GitBook Documentation Structure Tips](https://gitbook.com/docs/guides/docs-best-practices/documentation-structure-tips) - Best practices

### WebSearch (verified)
- [API Documentation Best Practices](https://konghq.com/blog/learning-center/guide-to-api-documentation) - Kong Inc.
- [Atlassian Documentation Best Practices](https://www.atlassian.com/blog/loom/software-documentation-best-practices)

## Metadata

**Confidence breakdown:**
- Code patterns: HIGH - Direct codebase analysis
- Documentation structure: HIGH - Industry-standard Diataxis framework
- Best practices: MEDIUM - WebSearch verified with multiple sources

**Research date:** 2026-02-01
**Valid until:** 2026-03-01 (stable domain, documentation practices evolve slowly)
