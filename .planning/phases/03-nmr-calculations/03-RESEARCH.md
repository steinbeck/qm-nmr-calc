# Phase 3: NMR Calculations - Research

**Researched:** 2026-01-19
**Domain:** NMR shielding calculation via ISiCLE/NWChem, conformer handling, presets
**Confidence:** HIGH

## Summary

Phase 3 extends the existing geometry optimization pipeline to produce NMR chemical shifts. Research confirms:

1. **ISiCLE provides a complete NMR pipeline.** The `isicle.qm.dft()` function accepts `tasks=['shielding']` to calculate NMR shielding constants via NWChem. The `NWChemParser` extracts isotropic shielding values with atom indices and atom types. Shielding-to-shift conversion requires linear regression parameters (slope/intercept) per nucleus type.

2. **Conformer handling is built into ISiCLE.** The `isicle.md.md()` function with CREST/xTB generates conformers, and `isicle.conformers.ConformationalEnsemble.reduce()` performs Boltzmann-weighted averaging of properties. However, conformer generation requires xTB/CREST binaries, which adds deployment complexity. A single-conformer option is simpler and often sufficient for rigid molecules.

3. **Two presets are straightforward to implement.** Draft uses smaller basis set (6-31G*) for speed; Production uses larger basis set (6-311+G(2d,p)) for accuracy. Both use B3LYP functional, which is well-validated for NMR. Solvation is mandatory via NWChem COSMO with user-specified solvent.

**Primary recommendation:** Extend `run_geometry_optimization()` to include NMR shielding calculation as a second DFT step. Store raw shielding values in JSON, convert to chemical shifts using published scaling factors. Implement presets as parameter dictionaries that configure basis set and task lists.

## Standard Stack

The established libraries/tools for NMR calculations:

### Core (Already Installed)

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| ISiCLE | 2.0.0 (fork) | NWChem wrapper, shielding extraction | Already handles NMR workflow natively |
| NWChem | 7.0.2 | DFT shielding calculation | GIAO NMR shielding via COSMO solvation |
| RDKit | >=2025.3.0 | Atom indexing, molecular properties | Map shielding values to atoms |

### Supporting (Already Installed)

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| numpy | (via ISiCLE) | Numerical operations | Boltzmann averaging |
| pandas | (via ISiCLE) | Data manipulation | Shielding data storage |

### Optional for Conformers

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| xTB | >=6.5 | Semi-empirical optimization | Conformer generation via CREST |
| CREST | >=3.0 | Conformer search | Automated conformer sampling |

**Note:** xTB/CREST are NOT currently installed. Single-conformer mode is the baseline; conformer averaging can be added later.

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| NWChem GIAO | ORCA NMR | ORCA not installed, ISiCLE supports both |
| B3LYP functional | mPW1PW91 | mPW1PW91 slightly better but B3LYP well-validated |
| COSMO solvation | SMD | SMD not available in NWChem for shielding |

## Architecture Patterns

### Recommended Extension to Project Structure

```
src/
└── qm_nmr_calc/
    ├── isicle_wrapper.py     # EXTEND: Add run_nmr_calculation()
    ├── tasks.py              # EXTEND: Add run_nmr_task()
    ├── models.py             # EXTEND: Add preset enum, NMR result models
    ├── storage.py            # EXTEND: Store NMR results
    ├── presets.py            # NEW: Preset configurations
    └── shifts.py             # NEW: Shielding-to-shift conversion
```

### Pattern 1: Two-Step DFT Workflow

**What:** Geometry optimization followed by NMR shielding calculation.

**When to use:** Always - NMR accuracy depends on optimized geometry.

**Example:**
```python
# Source: ISiCLE examples/nmr_chemical_shifts.py
import isicle

def run_nmr_calculation(smiles: str, job_dir: Path, preset: dict) -> dict:
    """Run geometry optimization + NMR shielding."""

    # 1. Load and embed molecule
    geom = isicle.load(smiles)
    geom = geom.initial_optimize(embed=True, forcefield='UFF', ff_iter=200)

    scratch_dir = job_dir / 'scratch'
    scratch_dir.mkdir(exist_ok=True)

    # 2. DFT geometry optimization
    opt_wrapper = isicle.qm.dft(
        geom,
        backend='NWChem',
        tasks=['optimize'],
        functional=preset['functional'],
        basis_set=preset['basis_set'],
        cosmo=True,
        solvent=preset['solvent'],
        scratch_dir=str(scratch_dir),
        processes=preset['processes'],
    )
    opt_result = opt_wrapper.parse()
    optimized_geom = opt_result['geometry']

    # 3. NMR shielding calculation on optimized geometry
    nmr_wrapper = isicle.qm.dft(
        optimized_geom,
        backend='NWChem',
        tasks=['shielding'],  # NMR shielding via GIAO
        functional=preset['functional'],
        basis_set=preset['nmr_basis_set'],  # Can be different for NMR
        cosmo=True,
        solvent=preset['solvent'],
        scratch_dir=str(scratch_dir),
        processes=preset['processes'],
    )
    nmr_result = nmr_wrapper.parse()

    return {
        'geometry': optimized_geom,
        'energy': nmr_result.get('energy'),
        'shielding': nmr_result.get('shielding'),
    }
```

### Pattern 2: Preset Configuration Dictionaries

**What:** Named presets that configure DFT parameters.

**When to use:** All job submissions - user selects preset name.

**Example:**
```python
# Source: User decisions in CONTEXT.md + DFT literature
from enum import Enum
from typing import TypedDict

class PresetName(str, Enum):
    DRAFT = "draft"
    PRODUCTION = "production"

class CalculationPreset(TypedDict):
    name: str
    description: str
    functional: str
    basis_set: str       # For geometry optimization
    nmr_basis_set: str   # For NMR shielding (can be same or different)
    processes: int
    max_iter: int

PRESETS: dict[PresetName, CalculationPreset] = {
    PresetName.DRAFT: {
        'name': 'draft',
        'description': 'Fast calculations for quick checks',
        'functional': 'b3lyp',
        'basis_set': '6-31G*',
        'nmr_basis_set': '6-31G*',
        'processes': 4,
        'max_iter': 100,
    },
    PresetName.PRODUCTION: {
        'name': 'production',
        'description': 'Balanced accuracy and compute time (default)',
        'functional': 'b3lyp',
        'basis_set': '6-31G*',           # Optimize with smaller basis
        'nmr_basis_set': '6-311+G(2d,p)', # NMR with larger basis for accuracy
        'processes': 4,
        'max_iter': 150,
    },
}

DEFAULT_PRESET = PresetName.PRODUCTION
```

### Pattern 3: Shielding-to-Shift Conversion

**What:** Convert calculated shielding constants to chemical shifts using linear regression.

**When to use:** After NMR calculation completes.

**Reference:** ISiCLE's `conformers.transform()` function uses: `shift = m * shielding + b`

**Example:**
```python
# Source: ISiCLE conformers.py + literature values
# Scaling factors for B3LYP/6-311+G(2d,p) in CDCl3 (common NMR solvent)
# Reference: Pierens, J. Comput. Chem. 2014, 35, 1388-1394

SCALING_FACTORS = {
    'H': {'m': -1.0, 'b': 31.8},   # TMS reference ~31.8 ppm shielding
    'C': {'m': -1.0, 'b': 182.5},  # TMS reference ~182.5 ppm shielding
}

def shielding_to_shift(shielding_data: dict, scaling: dict = SCALING_FACTORS) -> dict:
    """
    Convert shielding values to chemical shifts.

    Input shielding_data format (from NWChemParser):
        {'index': [1, 2, 3...], 'atom': ['H', 'C', 'H'...], 'shielding': [29.1, 156.2, ...]}

    Returns:
        {'index': [1, 2, 3...], 'atom': ['H', 'C'...], 'shift': [2.7, 26.3, ...]}
    """
    shifts = []
    for idx, atom, shield in zip(
        shielding_data['index'],
        shielding_data['atom'],
        shielding_data['shielding']
    ):
        if atom in scaling:
            m, b = scaling[atom]['m'], scaling[atom]['b']
            shift = m * shield + b
            shifts.append({
                'index': idx,
                'atom': atom,
                'shielding': shield,
                'shift': round(shift, 2),
            })

    # Separate by nucleus type
    h_shifts = [s for s in shifts if s['atom'] == 'H']
    c_shifts = [s for s in shifts if s['atom'] == 'C']

    return {
        '1H': sorted(h_shifts, key=lambda x: x['shift'], reverse=True),
        '13C': sorted(c_shifts, key=lambda x: x['shift'], reverse=True),
    }
```

### Pattern 4: Solvent Validation

**What:** Validate user-provided solvent against NWChem COSMO supported solvents.

**When to use:** At job submission time.

**Example:**
```python
# Source: NWChem documentation - Solvation-Models.html
# Common NMR solvents (subset of full 180+ solvent list)
SUPPORTED_SOLVENTS = {
    'chcl3': 'Chloroform (CDCl3)',
    'dmso': 'Dimethylsulfoxide (DMSO-d6)',
    'methanol': 'Methanol (CD3OD)',
    'acetone': 'Acetone (acetone-d6)',
    'benzene': 'Benzene (C6D6)',
    'h2o': 'Water (D2O)',
    'acetntrl': 'Acetonitrile (CD3CN)',
    'dcm': 'Dichloromethane (CD2Cl2)',
    'thf': 'Tetrahydrofuran (THF-d8)',
    'pyridine': 'Pyridine (pyridine-d5)',
    'toluene': 'Toluene (toluene-d8)',
}

def validate_solvent(solvent: str) -> str | None:
    """
    Validate solvent name. Returns normalized name or None if invalid.
    """
    normalized = solvent.lower().strip()
    if normalized in SUPPORTED_SOLVENTS:
        return normalized
    return None
```

### Anti-Patterns to Avoid

- **Skipping geometry optimization:** NMR shielding is highly sensitive to geometry. Always optimize first.
- **Using gas-phase for NMR:** Solvation effects are significant for NMR. Always use COSMO.
- **Reporting raw shielding values:** Users expect chemical shifts (ppm), not shielding constants.
- **Hardcoding TMS shielding:** Use proper linear regression scaling; TMS position varies with method.
- **Running NMR with wrong basis for nucleus:** 1H and 13C may benefit from different basis sets for accuracy vs speed.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| NWChem input for NMR | Manual .nw file | ISiCLE's _configure_shielding() | Handles GIAO, COSMO, all atoms |
| Shielding extraction | Regex parsing | NWChemParser._parse_shielding() | Tested, handles all formats |
| Atom type identification | Manual lookup | ISiCLE stores 'atom' in shielding dict | Already extracted from output |
| Conformer averaging | Manual Boltzmann | isicle.conformers.reduce() | Handles energy weighting correctly |
| Solvent parameters | Lookup table | NWChem built-in 'solvent X' keyword | 180+ solvents pre-parameterized |

**Key insight:** ISiCLE's `qm.dft()` with `tasks=['shielding']` handles NMR calculation. The parser extracts shielding with atom indices. Only need thin wrapper for preset application and shift conversion.

## Common Pitfalls

### Pitfall 1: Shielding Index Off-by-One

**What goes wrong:** Atom indices in NWChem output don't match RDKit atom indices.

**Why it happens:** NWChem uses 1-based indexing; RDKit uses 0-based indexing.

**How to avoid:**
- ISiCLE's parser returns NWChem indices (1-based)
- Subtract 1 when mapping to RDKit atoms
- Store both indices for clarity: `nwchem_idx` and `rdkit_idx`

**Warning signs:** Wrong atoms labeled with shifts; off-by-one patterns in results.

### Pitfall 2: Missing Atoms in Shielding Output

**What goes wrong:** Not all atoms appear in NMR results.

**Why it happens:** NWChem calculates shielding for ALL atoms by default, but output parsing may miss some. Also, ISiCLE filters to C and H by default in some workflows.

**How to avoid:**
- Check output atom count matches molecule
- Verify both 1H and 13C atoms are present
- Use `atoms=['C', 'H']` parameter explicitly if filtering needed

**Warning signs:** Fewer shifts than expected; missing carbons or hydrogens.

### Pitfall 3: Large Molecules Fail Silently

**What goes wrong:** NWChem runs out of memory without clear error.

**Why it happens:** 6-311+G(2d,p) basis requires significant memory. NWChem default memory settings too low.

**How to avoid:**
- Set explicit memory: `mem_global=4000, mem_heap=400, mem_stack=1200` (MB)
- For production preset on larger molecules (>30 atoms), may need more
- Consider draft preset for initial testing

**Warning signs:** NWChem exits without shielding data; ARMCI errors in logs.

### Pitfall 4: Wrong Scaling Factors for Method

**What goes wrong:** Calculated shifts systematically off by several ppm.

**Why it happens:** Scaling factors depend on functional, basis set, AND solvent. Using factors from different method gives poor results.

**How to avoid:**
- Match scaling factors to exact method used
- For B3LYP/6-311+G(2d,p) in CHCl3, use published values from Pierens et al.
- If using different method, recalibrate or use method-specific literature values

**Warning signs:** All shifts offset by constant amount; R-squared to experimental is poor.

### Pitfall 5: Conformer Generation Dependencies

**What goes wrong:** Conformer search fails because xTB/CREST not installed.

**Why it happens:** ISiCLE's `md.md()` with CREST requires xTB and CREST binaries in PATH.

**How to avoid:**
- Default to single-conformer mode (no xTB/CREST required)
- Document xTB/CREST as optional dependencies for conformer averaging
- Detect missing binaries at startup if conformer mode requested

**Warning signs:** FileNotFoundError for xtb or crest; conformer search hangs.

## Code Examples

Verified patterns from ISiCLE source code and NWChem documentation:

### Extended JobInput Model with Preset and Solvent

```python
# Source: CONTEXT.md decisions
from pydantic import BaseModel, Field
from typing import Optional, Literal

class JobInput(BaseModel):
    """Input parameters for a calculation job."""
    smiles: str
    name: Optional[str] = None
    preset: Literal["draft", "production"] = "production"
    solvent: str = Field(..., description="NWChem COSMO solvent name (required)")
```

### NMR Result Model

```python
# Source: ISiCLE parse.py NWChemParser output format
from pydantic import BaseModel
from typing import List, Optional

class AtomShift(BaseModel):
    """Chemical shift for a single atom."""
    index: int               # NWChem 1-based index
    atom: str                # Element symbol (H, C)
    shielding: float         # Raw isotropic shielding (ppm)
    shift: float             # Chemical shift (ppm, relative to TMS)

class NMRResults(BaseModel):
    """NMR calculation results."""
    h1_shifts: List[AtomShift]   # 1H chemical shifts
    c13_shifts: List[AtomShift]  # 13C chemical shifts

    # Calculation metadata
    functional: str
    basis_set: str
    solvent: str
```

### Extended JobStatus with NMR Results

```python
# Source: Extending existing models.py
class JobStatus(BaseModel):
    """Complete status of a calculation job."""
    job_id: str
    status: Literal["queued", "running", "complete", "failed"]
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None

    input: JobInput
    isicle_version: str
    nwchem_version: str

    # Results (populated on completion)
    nmr_results: Optional[NMRResults] = None
    optimized_geometry_file: Optional[str] = None

    # Error info
    error_message: Optional[str] = None
    error_traceback: Optional[str] = None
```

### Complete NMR Task Implementation

```python
# Source: Extending existing tasks.py
@huey.task()
def run_nmr_task(job_id: str) -> dict:
    """Execute NMR calculation for a queued job."""
    from .presets import PRESETS, PresetName
    from .shifts import shielding_to_shift

    job_status = load_job_status(job_id)
    if job_status is None:
        raise ValueError(f"Job {job_id} not found")

    preset_name = PresetName(job_status.input.preset)
    preset = PRESETS[preset_name]
    solvent = job_status.input.solvent
    smiles = job_status.input.smiles
    job_dir = get_job_dir(job_id)

    # Load molecule
    geom = isicle.load(smiles)
    geom = geom.initial_optimize(embed=True, forcefield='UFF', ff_iter=200)

    scratch_dir = job_dir / 'scratch'
    scratch_dir.mkdir(exist_ok=True)

    # Step 1: Geometry optimization
    opt_wrapper = isicle.qm.dft(
        geom,
        backend='NWChem',
        tasks=['optimize'],
        functional=preset['functional'],
        basis_set=preset['basis_set'],
        cosmo=True,
        solvent=solvent,
        max_iter=preset['max_iter'],
        scratch_dir=str(scratch_dir),
        processes=preset['processes'],
    )
    opt_result = opt_wrapper.parse()
    optimized_geom = opt_result['geometry']

    # Save optimized geometry
    output_dir = job_dir / 'output'
    output_dir.mkdir(exist_ok=True)
    geom_file = output_dir / 'optimized.xyz'
    isicle.save(str(geom_file), optimized_geom)

    # Step 2: NMR shielding calculation
    nmr_wrapper = isicle.qm.dft(
        optimized_geom,
        backend='NWChem',
        tasks=['shielding'],
        functional=preset['functional'],
        basis_set=preset['nmr_basis_set'],
        cosmo=True,
        solvent=solvent,
        scratch_dir=str(scratch_dir),
        processes=preset['processes'],
    )
    nmr_result = nmr_wrapper.parse()

    # Convert shielding to shifts
    shielding_data = nmr_result.get('shielding')
    if shielding_data is None:
        raise RuntimeError("NMR calculation did not produce shielding data")

    shifts = shielding_to_shift(shielding_data)

    # Save results
    results_file = output_dir / 'nmr_results.json'
    results_file.write_bytes(orjson.dumps({
        'shielding': shielding_data,
        'shifts': shifts,
        'functional': preset['functional'],
        'basis_set': preset['nmr_basis_set'],
        'solvent': solvent,
    }, option=orjson.OPT_INDENT_2))

    return {
        'success': True,
        'job_id': job_id,
        'output_files': {
            'geometry': str(geom_file),
            'nmr_results': str(results_file),
        }
    }
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| TMS-referenced shifts | Linear regression scaling | 2010s | Better accuracy across range |
| Gas-phase NMR | Implicit solvation (COSMO) | Standard practice | Matches experimental conditions |
| Single geometry | Conformer-averaged | Research-grade | Better for flexible molecules |
| Large basis everywhere | Split basis (geom/NMR) | Performance optimization | Same accuracy, less time |

**Current best practices:**
- B3LYP or mPW1PW91 functional for NMR (well-validated)
- 6-311+G(2d,p) or similar for NMR shielding calculation
- Smaller basis (6-31G*) acceptable for geometry optimization step
- Always include implicit solvation matching experimental conditions

## Open Questions

Things that couldn't be fully resolved:

1. **Optimal scaling factors for our method/solvent combination**
   - What we know: ISiCLE example uses m=-0.9921, b=32.2773 for 1H; m=-0.9816, b=180.4295 for 13C
   - What's unclear: These may be for specific functional/basis/solvent combination
   - Recommendation: Start with ISiCLE defaults, validate against known compounds, adjust if needed

2. **xTB/CREST installation for conformer averaging**
   - What we know: Not currently installed; ISiCLE supports it natively
   - What's unclear: Installation complexity, impact on deployment
   - Recommendation: Defer to Phase 3.1 or later; single-conformer is baseline

3. **Memory requirements for larger molecules**
   - What we know: NWChem defaults often insufficient for production basis set
   - What's unclear: Exact thresholds by atom count
   - Recommendation: Start with generous memory (4GB global), monitor and adjust

4. **Error recovery for partial NMR failures**
   - What we know: Geometry optimization may succeed but NMR fail
   - What's unclear: Whether to save partial results
   - Recommendation: Save geometry if optimization succeeds; mark job failed if NMR fails

## Sources

### Primary (HIGH confidence)

- **ISiCLE source code** - Direct analysis of local fork:
  - `isicle/qm.py` - NWChemWrapper._configure_shielding(), task handling
  - `isicle/parse.py` - NWChemParser._parse_shielding() output format
  - `isicle/conformers.py` - Boltzmann averaging, transform() for scaling
  - `isicle/md.py` - CREST/xTB conformer generation
  - `examples/nmr_chemical_shifts.py` - Complete workflow example

- **NWChem Official Documentation**:
  - [Solvation Models](https://nwchemgit.github.io/Solvation-Models.html) - COSMO solvent list (180+ solvents)
  - [Properties](https://nwchemgit.github.io/Properties.html) - NMR shielding calculation

- **xTB Documentation**:
  - [Implicit Solvation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html) - ALPB solvent list for conformer generation

### Secondary (MEDIUM confidence)

- [DFT NMR chemical shift calculations](https://www.nature.com/articles/s41598-022-22900-y) - B3LYP functional validation
- [1H/13C scaling factors](https://royalsocietypublishing.org/doi/10.1098/rsos.210954) - Geometry optimization impact
- [Pierens scaling factors](https://onlinelibrary.wiley.com/doi/10.1002/jcc.23638) - Published scaling factors for various methods/solvents

### Tertiary (LOW confidence)

- WebSearch results for NMR DFT best practices - Cross-verified with ISiCLE implementation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - ISiCLE source code directly analyzed, NWChem docs verified
- Architecture: HIGH - Based on ISiCLE's working nmr_chemical_shifts.py example
- Preset parameters: MEDIUM - Based on literature, may need tuning for specific use case
- Pitfalls: MEDIUM - Based on docs + ISiCLE issues, some edge cases unverified

**Research date:** 2026-01-19
**Valid until:** 2026-02-19 (stable libraries, 30 days)

---

*Phase: 03-nmr-calculations*
*Research completed: 2026-01-19*
