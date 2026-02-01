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

**Additional libraries:**
- [3Dmol.js](#3dmoljs) - Interactive 3D viewer with shift labels
- [SmilesDrawer](#smilesdrawer) - Real-time 2D structure preview
- [CREST/xTB](#crestxtb-optional) - Optional high-accuracy conformer generation

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

---

## Huey Task Queue

Huey provides async job processing with SQLite persistence, enabling long-running NWChem calculations to run in background workers.

**Official docs:** [Huey Documentation](https://huey.readthedocs.io/)

### Integration Points

| Module | Components | Purpose |
|--------|------------|---------|
| `queue.py` | `huey` instance, signal handlers | Queue configuration and lifecycle |
| `tasks.py` | `run_nmr_task()`, `run_ensemble_nmr_task()` | Task definitions |

### Why Huey over Celery

The project uses Huey instead of the more common Celery + Redis stack:

1. **Simpler deployment** - No Redis server required
2. **SQLite backend** - Uses `./data/huey.db` for durability (single file)
3. **Single-process consumer** - NWChem calculations are CPU-bound; multiple workers don't help
4. **Lightweight** - Minimal dependencies compared to Celery's ecosystem

### Queue Configuration

The queue is configured in `queue.py` with signal handlers for job lifecycle management.

```python
# Source: src/qm_nmr_calc/queue.py, lines 1-89
from huey import SqliteHuey, signals

# SQLite-backed queue with fsync for durability
huey = SqliteHuey(
    'qm-nmr-calc',
    filename='./data/huey.db',
    fsync=True
)

@huey.signal(signals.SIGNAL_EXECUTING)
def on_task_start(signal, task):
    """Update job status when task starts executing."""
    job_id = task.args[0]
    update_job_status(job_id, status='running', started_at=datetime.utcnow())

@huey.signal(signals.SIGNAL_COMPLETE)
def on_task_complete(signal, task):
    """Update job status and send notification on success."""
    job_id = task.args[0]
    update_job_status(job_id, status='complete', completed_at=datetime.utcnow())
    # Send email notification if configured
    job_status = load_job_status(job_id)
    if job_status and job_status.input.notification_email:
        send_job_notification_sync(to_email=..., job_id=job_id, status="complete")

@huey.signal(signals.SIGNAL_ERROR)
def on_task_error(signal, task, exc=None):
    """Update job status with error details on failure."""
    job_id = task.args[0]
    error_msg = str(exc) if exc else 'Unknown error'
    update_job_status(job_id, status='failed', error_message=error_msg)
```

**Key configuration:**
- `SqliteHuey` - SQLite-backed storage (no Redis required)
- `fsync=True` - Ensures durability across crashes
- `filename='./data/huey.db'` - Database location in data directory

**Signal handlers:**
| Signal | When Fired | Action |
|--------|------------|--------|
| `SIGNAL_EXECUTING` | Task starts | Set job status to `running`, record `started_at` |
| `SIGNAL_COMPLETE` | Task succeeds | Set job status to `complete`, send email notification |
| `SIGNAL_ERROR` | Task fails | Set job status to `failed`, record error message |
| `SIGNAL_INTERRUPTED` | Graceful shutdown | Set job status to `failed` with shutdown message |

### Task Definition

Tasks use the `@huey.task()` decorator and track progress through step callbacks.

```python
# Source: src/qm_nmr_calc/tasks.py, lines 82-144
from .queue import huey

@huey.task()
def run_nmr_task(job_id: str) -> dict:
    """Execute NMR calculation for a queued job."""
    job_status = load_job_status(job_id)
    preset = PRESETS[PresetName(job_status.input.preset)]

    # Step 1: Geometry optimization
    start_step(job_id, "geometry_optimization", "Optimizing geometry")

    # Callback to switch step tracking when optimization completes
    def on_opt_complete():
        start_step(job_id, "nmr_shielding", "Computing NMR shielding")

    result = run_calculation(
        smiles=job_status.input.smiles,
        job_dir=get_job_dir(job_id),
        preset=preset,
        solvent=job_status.input.solvent,
        on_optimization_complete=on_opt_complete,
    )

    # Step 3: Post-processing
    start_step(job_id, "post_processing", "Generating results")

    # Convert shielding to shifts, generate visualizations...
    return {'success': True, 'job_id': job_id}
```

**Task pattern:**
- `job_id` is always the first argument (used by signal handlers to identify the job)
- Step progress is tracked via `start_step()` and `complete_current_step()` functions
- The `status.json` file in the job directory is updated at each pipeline stage
- Exceptions propagate to Huey, triggering `SIGNAL_ERROR` for automatic status update

### Consumer Operations

Start the Huey consumer to process queued jobs:

```bash
# Start consumer (single worker, thread-based)
huey_consumer qm_nmr_calc.queue.huey -w 1 -k thread

# Useful flags:
#   -w N      Number of workers (default: 1)
#   -k thread Use threading instead of multiprocessing
#   -v        Verbose output
#   -l /path  Log file location
```

**Single worker recommendation:** NWChem calculations spawn MPI processes that consume all available CPU cores. Running multiple Huey workers would cause resource contention. Use `-w 1` for production.

**Startup order:**
1. Start the FastAPI application first (`uvicorn qm_nmr_calc.api.main:app`)
2. Then start the Huey consumer (`huey_consumer qm_nmr_calc.queue.huey`)

The API can accept and queue jobs even when the consumer is down; they'll process when the consumer restarts.

---

## 3Dmol.js

3Dmol.js provides interactive 3D molecular visualization with chemical shift labels overlaid on the structure.

**Official docs:** [3Dmol.js Documentation](https://3dmol.csb.pitt.edu/)

### Integration Points

| Template | Features | Purpose |
|----------|----------|---------|
| `results.html` | Viewer, shift labels, conformer switching | Results visualization |
| `status.html` | Viewer (optional) | Progress preview |

**CDN:** `https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.3/3Dmol-min.js`

### Viewer Setup

The viewer is initialized in the results page to display the optimized 3D geometry.

```javascript
// Source: src/qm_nmr_calc/api/templates/results.html, lines 401-409
function initViewer() {
    const container = document.getElementById('viewer-container-3d');
    if (!container) return;

    viewer = $3Dmol.createViewer(container, {
        backgroundColor: 'white'
    });
    loadGeometryData();
}
```

**Key configuration:**
- Container element via `getElementById` (not jQuery selector)
- `backgroundColor: 'white'` for clean presentation
- Viewer instance stored globally for label/conformer operations

### Geometry Loading

Geometry data is fetched from the API and loaded with format preference.

```javascript
// Source: src/qm_nmr_calc/api/templates/results.html, lines 411-439
async function loadGeometryData() {
    const response = await fetch('/api/v1/jobs/' + JOB_ID + '/geometry.json');
    geometryData = await response.json();

    // ... conformer handling ...
    displayConformer(0);
}

function displayConformer(index) {
    // Clear and redraw viewer
    viewer.removeAllModels();
    viewer.removeAllLabels();

    // Add model (prefer SDF for proper bonds)
    let model;
    if (sdf) {
        model = viewer.addModel(sdf, 'sdf');
    } else if (xyz) {
        model = viewer.addModel(xyz, 'xyz', {assignBonds: true});
    }

    // Apply style
    viewer.setStyle({}, {
        stick: {radius: 0.12, colorscheme: 'Jmol'},
        sphere: {scale: 0.25, colorscheme: 'Jmol'}
    });

    viewer.zoomTo();
    viewer.render();
}
```

**Format preference:**
- **SDF preferred** - Contains explicit bond orders for correct aromatic ring and double bond visualization
- **XYZ fallback** - Uses `assignBonds: true` to infer bonds from distances (less accurate for aromatics)

**Style configuration:**
- `stick` - Bond representation with 0.12 radius
- `sphere` - Atom centers with 0.25 scale
- `colorscheme: 'Jmol'` - Standard Jmol element colors (carbon gray, oxygen red, etc.)

### Chemical Shift Labels

The `addShiftLabels()` function overlays NMR shift values directly on atoms in the 3D view.

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

**CRITICAL - Index Conversion:**
- 3Dmol.js `atom.serial` is **0-based** (first atom is 0)
- NWChem atom indices are **1-based** (first atom is 1)
- Always add 1 when looking up assignments: `atom.serial + 1`

**Label styling:**
| Element | Background Color | Hex Code |
|---------|-----------------|----------|
| 1H | Blue | `#3498db` |
| 13C | Orange | `#e67e22` |

**Label options:**
- `fontColor: 'white'` - High contrast on colored backgrounds
- `backgroundOpacity: 0.85` - Slightly transparent for depth perception
- `inFront: true` - Labels always visible (not occluded by molecule)
- `fontSize: 11` - Readable without obscuring structure

**Note:** Labels are added per-atom, not per-peak. For equivalent atoms with identical shifts (e.g., methyl group hydrogens), each atom gets its own label with the same value.

---

## SmilesDrawer

SmilesDrawer provides real-time 2D structure preview as users type SMILES strings in the submission form.

**Official docs:** [SmilesDrawer GitHub](https://github.com/reymond-group/smilesDrawer)

### Integration Points

| Template | Features | Purpose |
|----------|----------|---------|
| `submit.html` | Canvas preview, validation feedback | Input assistance |

**CDN:** `https://unpkg.com/smiles-drawer@1.2.0/dist/smiles-drawer.min.js`

### Canvas Setup

The drawer is configured with custom colors matching the application theme.

```javascript
// Source: src/qm_nmr_calc/api/templates/submit.html, lines 149-173
const smilesDrawer = new SmilesDrawer.Drawer({
    width: 350,
    height: 350,
    bondThickness: 0.6,
    bondLength: 15,
    terminalCarbons: true,
    explicitHydrogens: false,
    themes: {
        light: {
            C: '#222',
            O: '#e74c3c',
            N: '#3498db',
            S: '#f1c40f',
            P: '#e67e22',
            F: '#1abc9c',
            CL: '#1abc9c',
            BR: '#e74c3c',
            I: '#9b59b6',
            H: '#aaa'
        }
    }
});
```

**Key options:**
- `width/height: 350` - Canvas dimensions in pixels
- `bondThickness: 0.6` - Bond line width
- `terminalCarbons: true` - Show CH3 labels explicitly
- `explicitHydrogens: false` - Hide non-terminal hydrogens for clarity
- Custom theme colors matching the application color palette

### Debounced Input Preview

SMILES parsing runs on input change with debouncing to avoid excessive redraws.

```javascript
// Source: src/qm_nmr_calc/api/templates/submit.html, lines 191-214
function updateMoleculePreview(smiles) {
    setupCanvas();

    if (!smiles || smiles.trim() === '') {
        const ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        status.textContent = 'Enter a SMILES string to preview';
        status.className = 'preview-card__status preview-card__status--empty';
        return;
    }

    SmilesDrawer.parse(smiles, function(tree) {
        // Valid SMILES - draw molecule
        smilesDrawer.draw(tree, 'molecule-preview', 'light', false);
        status.textContent = 'Valid SMILES';
        status.className = 'preview-card__status preview-card__status--valid';
    }, function(err) {
        // Invalid SMILES - show error
        const ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        status.textContent = 'Invalid SMILES syntax';
        status.className = 'preview-card__status preview-card__status--invalid';
    });
}

// Debounced preview update (400ms delay)
const debouncedPreview = debounce(updateMoleculePreview, 400);
```

**Parse + Draw pattern:**
1. `SmilesDrawer.parse(smiles, successCallback, errorCallback)` - Parse SMILES string
2. On success: `smilesDrawer.draw(tree, canvasId, theme, isAnimated)` - Draw to canvas
3. On error: Clear canvas and show error status

**Status feedback:**
| Status | CSS Class | Message |
|--------|-----------|---------|
| Empty | `--empty` | "Enter a SMILES string to preview" |
| Valid | `--valid` | "Valid SMILES" |
| Invalid | `--invalid` | "Invalid SMILES syntax" |

**Debounce timing:** 400ms delay prevents rapid redraws during typing while remaining responsive.

---

## CREST/xTB (Optional)

CREST (Conformer-Rotamer Ensemble Sampling Tool) with xTB provides high-accuracy conformer ensemble generation as an optional alternative to RDKit.

**Official docs:**
- [CREST Documentation](https://crest-lab.github.io/crest-docs/)
- [xTB Documentation](https://xtb-docs.readthedocs.io/)

### Integration Points

| Module | Functions | Purpose |
|--------|-----------|---------|
| `conformers/crest_generator.py` | `detect_crest_available()`, `run_crest()`, `parse_crest_ensemble()` | Ensemble generation |

**Binaries required:** Both `crest` and `xtb` must be available on PATH.

**Fallback behavior:** When CREST is not available, the system automatically falls back to RDKit KDG conformer generation.

### Availability Detection

CREST availability is checked once at startup and cached.

```python
# Source: src/qm_nmr_calc/conformers/crest_generator.py, lines 47-64
import shutil
from functools import lru_cache

@lru_cache(maxsize=1)
def detect_crest_available() -> bool:
    """
    Detect if CREST and xTB binaries are available on PATH.

    CREST requires both binaries to function. This function uses a cache
    since binary availability doesn't change during process lifetime.

    Returns:
        True if both crest and xtb are found, False otherwise
    """
    crest_found = shutil.which("crest") is not None
    xtb_found = shutil.which("xtb") is not None
    return crest_found and xtb_found
```

**Pattern notes:**
- `@lru_cache(maxsize=1)` - Cache result (binaries don't appear/disappear during runtime)
- `shutil.which()` - Cross-platform binary lookup (like Unix `which` command)
- Both binaries required - CREST orchestrates, xTB performs calculations

**Health endpoint:** The `/api/v1/health` endpoint reports CREST availability in its response:
```json
{"crest_available": true}
```

### CREST Execution

CREST is run as a subprocess with timeout handling for complex molecules.

```python
# Source: src/qm_nmr_calc/conformers/crest_generator.py, lines 172-259
def run_crest(
    input_xyz: Path,
    solvent: str,
    charge: int = 0,
    ewin_kcal: float = 6.0,
    num_threads: int | None = None,
    timeout_seconds: int = 7200,
) -> Path:
    """Run CREST conformational search subprocess."""
    # Build CREST command
    cmd = [
        "crest",
        str(input_xyz),
        "--gfn2",           # GFN2-xTB method
        "--alpb", solvent,  # ALPB implicit solvent
        "--chrg", str(charge),
        "--ewin", str(ewin_kcal),
        "-T", str(num_threads),
    ]

    # Set up environment variables for stability
    env = os.environ.copy()
    env["OMP_STACKSIZE"] = "2G"
    env["GFORTRAN_UNBUFFERED_ALL"] = "1"

    subprocess.run(
        cmd,
        cwd=input_xyz.parent,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        env=env,
        check=True,
    )

    return input_xyz.parent / "crest_conformers.xyz"
```

**Key parameters:**
- `--gfn2` - Use GFN2-xTB method (good accuracy/speed balance)
- `--alpb` - ALPB implicit solvation model
- `--ewin` - Energy window in kcal/mol for conformer inclusion
- `-T` - Number of threads (defaults to all available)

**Environment variables:**
- `OMP_STACKSIZE=2G` - Increase OpenMP stack for large molecules
- `GFORTRAN_UNBUFFERED_ALL=1` - Prevent output buffering issues

**Timeout handling:**
- Default: 7200 seconds (2 hours)
- Macrocycles and complex molecules may require longer
- On timeout: `RuntimeError` raised (no silent fallback to RDKit)

**Solvent mapping (ALPB):**
| Job Solvent | ALPB Identifier |
|-------------|-----------------|
| `chcl3` | `chcl3` |
| `dmso` | `dmso` |
| `vacuum` | (no solvent flag) |

### Ensemble Parsing

CREST outputs conformers as concatenated XYZ files with energies in the comment line.

```python
# Source: src/qm_nmr_calc/conformers/crest_generator.py, lines 98-169
def parse_crest_ensemble(ensemble_file: Path) -> list[CRESTConformer]:
    """
    Parse CREST multi-structure XYZ ensemble file.

    CREST outputs conformers in concatenated XYZ format:
    - Each structure starts with atom count line
    - Comment line contains energy as first token (in Hartree)
    - Coordinate lines follow
    - Structures are concatenated sequentially
    """
    conformers = []
    lines = ensemble_file.read_text().strip().split("\n")

    i = 0
    conf_number = 1

    while i < len(lines):
        num_atoms = int(lines[i].strip())
        i += 1
        comment_line = lines[i].strip()
        energy_hartree = float(comment_line.split()[0])
        i += 1

        coord_lines = []
        for _ in range(num_atoms):
            coord_lines.append(lines[i])
            i += 1

        # Build xyz_block
        xyz_block = f"{num_atoms}\n{comment_line}\n"
        xyz_block += "\n".join(coord_lines) + "\n"

        conformer = CRESTConformer(
            conformer_id=f"conf_{conf_number:03d}",
            energy_hartree=energy_hartree,
            xyz_block=xyz_block,
        )
        conformers.append(conformer)
        conf_number += 1

    return conformers
```

**CREST output format:**
```
12                      <- atom count
-12.345678901234        <- energy (Hartree) in comment line
C   0.000000   0.000000   0.000000
H   1.089000   0.000000   0.000000
...
12                      <- next conformer
-12.345612345678
...
```

**CRESTConformer dataclass:**
| Field | Type | Description |
|-------|------|-------------|
| `conformer_id` | `str` | Sequential ID (e.g., `"conf_001"`) |
| `energy_hartree` | `float` | Total energy in Hartree |
| `xyz_block` | `str` | Complete XYZ content for this conformer |

**Energy handling:**
- CREST outputs absolute energies in Hartree
- After parsing, energies are converted to relative kcal/mol for filtering and Boltzmann weighting
- Conversion factor: 1 Hartree = 627.509 kcal/mol

---

**Navigation:** [Documentation Index](README.md) | [Technical Architecture](architecture.md)

**Dependencies:** See [pyproject.toml](../pyproject.toml) for library versions.

