# Pitfalls Research: Conformational Sampling for NMR Predictions

**Domain:** Conformational Sampling for NMR Predictions
**Researched:** 2026-01-26
**Confidence:** HIGH

---

## Critical Pitfalls

### Pitfall 1: Boltzmann Weight Numerical Overflow/Underflow

**What goes wrong:**
When calculating `exp(-ΔE/RT)` for conformers with large energy differences (>10 kcal/mol), numerical overflow or underflow occurs. With energy differences of 15+ kcal/mol, `exp(-15/0.59) ≈ 1e-11` causes float underflow to exactly 0.0, breaking normalization (division by zero). Conversely, programming errors like `exp(+ΔE/RT)` (wrong sign) cause instant overflow for any positive energy.

**Why it happens:**
- Naive implementation of Boltzmann formula without numerical stabilization
- Large conformer ensembles span wide energy ranges (0-20 kcal/mol common with CREST)
- Python's `math.exp()` has limited range: overflow at ~710, underflow at ~-745

**How to avoid:**
Use the **log-sum-exp trick** (exp-normalize):
```python
import numpy as np

def boltzmann_weights(energies_kcal, temperature=298.15):
    """Numerically stable Boltzmann weighting."""
    R = 0.001987204  # kcal/(mol·K)
    RT = R * temperature

    # Shift by minimum energy (makes lowest energy = 0)
    rel_energies = energies_kcal - np.min(energies_kcal)

    # Exponentiate shifted values (max exponent is 0, no overflow)
    boltz_factors = np.exp(-rel_energies / RT)

    # Normalize
    return boltz_factors / np.sum(boltz_factors)
```

**Warning signs:**
- All Boltzmann weights become 0.0 except one (underflow)
- `ZeroDivisionError` in normalization
- `RuntimeWarning: overflow encountered in exp`
- Weighted average exactly equals single conformer shift

**Phase to address:** Phase 2 (Boltzmann averaging implementation)

**Recovery strategy:**
If already deployed with naive implementation, filter conformers to energy window (ΔE < 10 kcal/mol) before Boltzmann weighting as temporary mitigation. Add numerical tests with wide energy ranges (0-20 kcal/mol).

---

### Pitfall 2: Wrong Energy Units for Boltzmann Weighting

**What goes wrong:**
Mixing energy units between conformer generation (xTB outputs Hartrees), DFT optimization (NWChem outputs Hartrees), and Boltzmann formula (expects kcal/mol or kJ/mol). Using Hartrees directly in `exp(-ΔE/RT)` with `RT = 0.59 kcal/mol` gives nonsensical weights (off by factor of 627).

**Why it happens:**
- CREST `crest_conformers.xyz` contains energies in **Hartrees** in comment lines
- NWChem DFT optimization outputs energies in **Hartrees** (search for "Total DFT energy")
- Boltzmann constant `R = 0.001987 kcal/(mol·K)` expects energies in **kcal/mol**
- Easy to forget conversion when parsing NWChem output

**How to avoid:**
Enforce unit consistency at parsing stage:
```python
# Constants
HARTREE_TO_KCAL = 627.5095  # CODATA 2022 value

def parse_nwchem_energy(output_text: str) -> float:
    """Extract DFT energy in kcal/mol from NWChem output."""
    # Find "Total DFT energy = -156.123456789" (Hartree)
    match = re.search(r'Total DFT energy\s*=\s*([-\d.]+)', output_text)
    if not match:
        raise ValueError("Could not find Total DFT energy in NWChem output")

    energy_hartree = float(match.group(1))
    return energy_hartree * HARTREE_TO_KCAL  # Convert to kcal/mol
```

Add explicit unit tracking in conformer data structures:
```python
@dataclass
class ConformerData:
    geometry_xyz: str
    energy_kcal: float  # Always kcal/mol after parsing
    source: str  # 'crest', 'rdkit_mmff', 'nwchem_dft'
```

**Warning signs:**
- All conformers have essentially identical Boltzmann weights (energy differences too small)
- Boltzmann weights don't match expected Boltzmann distribution (e.g., 3 kcal/mol difference should give ~100:1 ratio)
- Averaged NMR shifts are nonsensical (all conformers weighted equally despite large energy differences)

**Phase to address:** Phase 3 (DFT energy extraction)

**Recovery strategy:**
Add unit tests with known energy differences. Example: Two conformers at 0 and 1 kcal/mol should give weights ~0.843 and ~0.157 at 298K. If test fails, check unit conversion.

---

### Pitfall 3: CREST Hanging/Timeout on Certain Molecules

**What goes wrong:**
CREST enters infinite loop or runs for days on specific molecular structures, particularly macrocycles (>12-membered rings), molecules with many rotatable bonds (>10), or highly symmetric systems. No progress output, process becomes unresponsive, eventually hits timeout or fills disk with metadynamics trajectories.

**Why it happens:**
- CREST uses metadynamics for conformer sampling, which can get trapped in high-energy regions for complex molecules
- Initial geometry optimization failure causes CREST to retry indefinitely ([GitHub Issue #105](https://github.com/crest-lab/crest/issues/105), [Issue #313](https://github.com/crest-lab/crest/issues/313))
- Known issue on macOS 12 with GCC compiler ([Issue #382](https://github.com/crest-lab/crest/issues/382))
- Crude pre-optimization step can have 0% success rate on certain structures ([Issue #242](https://github.com/crest-lab/crest/issues/242))

**How to avoid:**
Implement multi-layer timeout and progress monitoring:
```python
import subprocess
from pathlib import Path

def run_crest_with_timeout(
    xyz_file: Path,
    timeout_seconds: int = 3600,  # 1 hour default
    progress_check_interval: int = 300,  # 5 minutes
):
    """Run CREST with timeout and progress monitoring."""
    cmd = ['crest', str(xyz_file), '--gfn2', '--T', '4']

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
        if process.returncode != 0:
            raise RuntimeError(f"CREST failed: {stderr}")
        return parse_crest_output('crest_conformers.xyz')

    except subprocess.TimeoutExpired:
        process.kill()
        raise TimeoutError(
            f"CREST exceeded {timeout_seconds}s timeout. "
            "Consider using RDKit instead for this molecule."
        )
```

Set environment variable for macOS compatibility:
```python
import os
os.environ['GFORTRAN_UNBUFFERED_ALL'] = '1'  # Fix macOS 12 issues
```

**Warning signs:**
- CREST process running >1 hour with no new conformers generated
- `/tmp` or scratch directory filling with large trajectory files
- No console output for >10 minutes
- CPU usage drops to 0% but process still running (hung state)

**Phase to address:** Phase 4 (CREST integration)

**Recovery strategy:**
If CREST times out, fall back to RDKit conformer generation automatically:
```python
try:
    conformers = generate_conformers_crest(mol, timeout=3600)
except TimeoutError:
    logger.warning("CREST timeout, falling back to RDKit")
    conformers = generate_conformers_rdkit(mol, num_confs=50)
```

---

### Pitfall 4: Atom Ordering Inconsistency Across Conformers

**What goes wrong:**
Atom indices change between conformer generation, DFT optimization, and NMR calculation. CREST outputs XYZ with atoms sorted by element type. NWChem might reorder atoms during geometry optimization. When averaging NMR shifts, atom `index=5` in conformer 1 might be carbon, but hydrogen in conformer 2, producing nonsensical averaged shifts.

**Why it happens:**
- Different programs use different internal atom ordering conventions
- CREST cregen deduplication can reorder atoms ([GitHub Issue #242](https://github.com/crest-lab/crest/issues/242))
- NWChem's `noautoz noautosym` prevents Z-matrix conversion but doesn't guarantee atom order preservation during optimization
- When converting between formats (SMILES → RDKit mol → XYZ → NWChem → XYZ), atom order can change at each step
- Recent research found "9 out of 12 structures examined had significant structural differences, including different orientations of functional groups" when comparing chiral molecule conformers ([Nature Communications 2025](https://www.nature.com/articles/s41467-025-66247-0))

**How to avoid:**
Establish canonical atom ordering at the start and enforce through all steps:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def canonical_atom_order_smiles(smiles: str) -> str:
    """Generate canonical SMILES with consistent atom ordering."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    # RDKit's canonical SMILES ensures reproducible atom order
    return Chem.MolToSmiles(mol, canonical=True, allHsExplicit=True)

def verify_atom_consistency(conformers: list[dict]) -> None:
    """Verify all conformers have identical atom sequences."""
    reference_atoms = conformers[0]['atom_sequence']  # ['C', 'C', 'H', 'H', ...]

    for i, conf in enumerate(conformers[1:], start=1):
        if conf['atom_sequence'] != reference_atoms:
            raise ValueError(
                f"Conformer {i} has different atom ordering than reference. "
                f"Expected: {reference_atoms[:5]}..., "
                f"Got: {conf['atom_sequence'][:5]}..."
            )
```

Track atom identity through all steps:
```python
@dataclass
class AtomIdentifier:
    """Immutable atom identity preserved across conformers."""
    element: str  # 'C', 'H', etc.
    canonical_index: int  # Index from canonical SMILES

@dataclass
class ConformerGeometry:
    atoms: list[AtomIdentifier]  # Fixed across all conformers
    coords: np.ndarray  # (N_atoms, 3), varies per conformer
```

**Warning signs:**
- Averaged 1H shift is in 13C range (100+ ppm) or vice versa
- Number of atoms differs between conformers
- Averaged shift for "atom 5" has huge variance (>50 ppm) across conformers
- Element type changes for same atom index

**Phase to address:** Phase 1 (Conformer data model design) and Phase 5 (NMR averaging)

**Recovery strategy:**
If inconsistency detected post-deployment:
1. Parse element type from NWChem output for each conformer
2. Re-index atoms by canonical SMILES ordering
3. Group shifts by (canonical_index, element) before averaging

---

### Pitfall 5: NWChem Scratch Directory Conflicts in Parallel Conformer Calculations

**What goes wrong:**
Running multiple NWChem jobs in parallel for different conformers causes file collisions in scratch directory. NWChem uses fixed filenames like `molecule.db`, `molecule.movecs`, `molecule.gridpts`. Multiple processes overwrite each other's files, causing "database corrupted" errors, convergence failures, or returning wrong results for wrong conformers.

**Why it happens:**
- Current code uses shared scratch directory: `job_dir / "scratch"` ([runner.py:138](https://github.com))
- NWChem's `start molecule` directive creates files prefixed with "molecule"
- Single-conformer code reuses same scratch directory, which works fine
- Parallel conformer calculations need isolated scratch spaces
- NWChem scratch directory guidance suggests "different jobs can share the same scratch directory as long as they use different prefix names" ([NWChem Wiki](https://github.com/nwchemgit/nwchem/wiki/Scratch_Dir))

**How to avoid:**
Create unique scratch directory per conformer using conformer ID:

```python
def run_calculation_conformer(
    smiles: str,
    job_dir: Path,
    conformer_id: int,  # Add this parameter
    preset: dict,
    solvent: str,
    processes: int = 4,
) -> dict:
    """Run NMR calculation for a single conformer."""
    # Unique scratch directory per conformer
    scratch_dir = job_dir / "scratch" / f"conf_{conformer_id}"
    scratch_dir.mkdir(parents=True, exist_ok=True)

    output_dir = job_dir / "output" / f"conf_{conformer_id}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Unique NWChem start name prevents database conflicts
    nwchem_start_name = f"conf_{conformer_id}"

    # Modify input generation to use unique start name:
    # start conf_0  (instead of start molecule)
    opt_input = generate_optimization_input(
        geometry_xyz=xyz_block,
        functional=preset["functional"],
        basis_set=preset["basis_set"],
        solvent=solvent,
        max_iter=preset.get("max_iter", 150),
        start_name=nwchem_start_name,  # NEW
    )

    # Rest of calculation...
```

**Warning signs:**
- "Database file is corrupted" errors in NWChem output
- Intermittent convergence failures (works sometimes, fails other times)
- Optimized geometry for conformer N looks identical to conformer M
- File `molecule.db` has modification time during conformer calculation
- Different shielding values when re-running same conformer

**Phase to address:** Phase 3 (Multi-conformer NWChem integration)

**Recovery strategy:**
If already deployed without isolation:
1. CRITICAL: Never run conformers in parallel (sequential only)
2. Add cleanup step between conformers: `rm -f scratch/molecule.*`
3. Longer-term: Refactor to per-conformer scratch directories

---

### Pitfall 6: Using Wrong Energy Level for Boltzmann Weighting

**What goes wrong:**
Weighting conformers by MMFF94 energies from RDKit, or GFN2-xTB energies from CREST, instead of DFT energies from NWChem optimization. Force field energies have poor correlation with DFT (Spearman ρ ~ -0.1 to -0.45), causing wrong conformer populations and wrong averaged NMR shifts. A conformer that's high-energy at MMFF level might be low-energy at DFT level.

**Why it happens:**
- RDKit `MMFFOptimizeMoleculeConfs()` returns MMFF energies, tempting to use immediately
- CREST outputs GFN2-xTB energies in `crest_conformers.xyz`, already computed
- DFT energies require parsing NWChem output after each optimization
- Reference analysis shows "MMFF energy ranking has poor correlation with DFT (Spearman ρ ~ -0.1 to -0.45)" ([conformational_sampling_nmr_analysis.md](https://github.com))
- GFN2-xTB is better (ρ ~ 0.39-0.47) but still not DFT-accurate

**How to avoid:**
Always use DFT energies from NWChem geometry optimization for Boltzmann weighting:

```python
@dataclass
class ConformerData:
    geometry_xyz: str
    mmff_energy_kcal: Optional[float] = None  # For filtering only
    xtb_energy_kcal: Optional[float] = None   # For filtering only
    dft_energy_kcal: Optional[float] = None   # For Boltzmann weighting

def boltzmann_average_shifts(conformers: list[ConformerData]) -> dict:
    """Average NMR shifts using DFT energies only."""
    # Verify all conformers have DFT energies
    if any(c.dft_energy_kcal is None for c in conformers):
        raise ValueError(
            "All conformers must have DFT energies for Boltzmann weighting. "
            "Run NWChem geometry optimization on all conformers first."
        )

    dft_energies = np.array([c.dft_energy_kcal for c in conformers])
    weights = boltzmann_weights(dft_energies, temperature=298.15)

    # ... rest of averaging
```

Use MMFF/xTB energies only for early filtering:
```python
def filter_conformers_by_mmff(conformers: list, window_kcal: float = 10.0):
    """Pre-filter conformers by MMFF energy before expensive DFT.

    NOTE: This is a rough filter. Final Boltzmann weights use DFT energies.
    """
    min_energy = min(c.mmff_energy_kcal for c in conformers)
    return [c for c in conformers if (c.mmff_energy_kcal - min_energy) <= window_kcal]
```

**Warning signs:**
- Averaged NMR shift deviates significantly from experimental (>5 ppm for 1H)
- Minor conformer by MMFF becomes major conformer by DFT
- Re-ranking conformers by DFT energy shows different order than MMFF/xTB
- Code calculates Boltzmann weights before DFT optimization completes

**Phase to address:** Phase 3 (DFT optimization) and Phase 5 (Boltzmann averaging)

**Recovery strategy:**
If deployed with wrong energy source:
1. Re-extract DFT energies from existing NWChem output files
2. Recalculate Boltzmann weights using DFT energies
3. Re-average NMR shifts with corrected weights
4. Add validation: compare MMFF vs DFT conformer rankings, log warning if very different

---

## Technical Debt Patterns

### Pattern 1: Insufficient Conformer Sampling

**What it looks like:**
Using fixed `num_confs=20` for all molecules, regardless of flexibility. Rigid molecules (0-1 rotatable bonds) get oversampled (wasting compute). Flexible molecules (8+ rotatable bonds) get undersampled (missing important conformers, poor NMR prediction accuracy).

**Why it happens:**
- Easy to hard-code a single conformer count
- No obvious failure signal (code runs, produces results)
- Only detectable by comparing to experimental NMR (requires validation dataset)

**Prevention:**
Adaptive conformer count based on molecular flexibility:
```python
def calculate_conformer_count(mol: Chem.Mol) -> int:
    """Calculate appropriate conformer count based on flexibility."""
    n_rotatable = Chem.Descriptors.NumRotatableBonds(mol)

    if n_rotatable <= 1:
        return 5  # Rigid, few conformers needed
    elif n_rotatable <= 4:
        return 20  # Semi-flexible
    elif n_rotatable <= 8:
        return 50  # Flexible
    else:
        return 100  # Very flexible, needs extensive sampling
```

**Detection:**
- Compare predicted vs experimental NMR shifts on benchmark set
- Large deviations (MAE > 1 ppm for 1H) suggest sampling issues
- Check conformer RMSD matrix: if all conformers are similar (RMSD < 1 Å), undersampled

---

### Pattern 2: No Conformer Deduplication

**What it looks like:**
RDKit generates 100 conformers, but 30 are duplicates (RMSD < 0.5 Å). Running DFT on all 100 wastes 30% of compute time. Boltzmann averaging over-weights the duplicated conformer structure (degeneracy issue).

**Why it happens:**
- RDKit's `pruneRmsThresh` parameter not set, or set too low
- Generated conformers not explicitly deduplicated after generation
- Default `pruneRmsThresh` might miss near-duplicates

**Prevention:**
Explicit deduplication after conformer generation:
```python
from rdkit.Chem import AllChem, rdMolAlign

def deduplicate_conformers(mol: Chem.Mol, rms_threshold: float = 0.5) -> list[int]:
    """Deduplicate conformers by RMSD, return unique conformer IDs."""
    conformer_ids = list(range(mol.GetNumConformers()))
    unique_ids = []

    for cid in conformer_ids:
        is_duplicate = False
        for uid in unique_ids:
            rms = rdMolAlign.AlignMol(mol, mol, prbCid=cid, refCid=uid)
            if rms < rms_threshold:
                is_duplicate = True
                break

        if not is_duplicate:
            unique_ids.append(cid)

    return unique_ids
```

**Detection:**
- Log conformer count before and after deduplication
- Calculate pairwise RMSD matrix, look for off-diagonal values < 0.5 Å
- Monitor DFT optimization: if many conformers converge to identical geometries, duplication occurred

---

### Pattern 3: RDKit ETKDG vs KDG Selection

**What it looks like:**
Using `ETKDGv3` for solution-phase NMR predictions. ETKDG biases toward crystal structure torsion angles (from CSD database), which differs from solution-phase conformer distributions. Predictions systematically deviate from experimental NMR.

**Why it happens:**
- ETKDG is the default in many tutorials and examples
- Not obvious that ETKDG has crystal structure bias
- Reference states "ETKDGv3 uses CSD-derived torsion preferences (biased toward crystal structures). For solution-phase NMR, pure distance geometry (KDG) may sample better." ([conformational_sampling_nmr_analysis.md](https://github.com))

**Prevention:**
Use KDG (Knowledge-free Distance Geometry) for solution-phase NMR:
```python
def generate_conformers_rdkit(mol: Chem.Mol, num_confs: int = 50):
    """Generate conformers for solution-phase NMR prediction."""
    mol = Chem.AddHs(mol)

    # KDG for solution sampling (no crystal bias)
    params = AllChem.KDGParams()
    params.pruneRmsThresh = 0.5
    params.numThreads = 0

    # NOT ETKDGv3() - that's for crystal/docking applications
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    # MMFF optimization
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)

    return mol
```

**Detection:**
- Compare predictions to experimental NMR for solution-phase measurements
- If ETKDG predictions systematically better for solid-state NMR than solution NMR, confirms crystal bias

---

## Integration Gotchas

### Gotcha 1: xTB Memory Explosion on Large Molecules

**Problem:**
xTB killed by OOM (out-of-memory) when CREST runs on molecules >500 atoms, or proteins. RAM usage spikes to 16+ GB, process terminated. Search results show "xTB gets killed when RAM usage goes to about 16 GB" for 5-40 kDa proteins ([GitHub Discussion #811](https://github.com/grimme-lab/xtb/discussions/811)).

**Root cause:**
- xTB allocates large matrices for semi-empirical calculations
- Stack overflow: "Some parts of the xtb program can be quite wasteful with stack memory" ([GitHub Issue #700](https://github.com/grimme-lab/xtb/issues/700))
- Systems >831-833 atoms always crash ([GitHub Issue #439](https://github.com/grimme-lab/xtb/issues/439))
- WSL2 limits RAM to 16 GB by default

**Mitigation:**
Set environment variables before running CREST:
```python
import os

def setup_xtb_environment():
    """Configure environment for xTB/CREST to avoid memory issues."""
    # Increase stack size (Linux/macOS)
    os.environ['OMP_STACKSIZE'] = '15G'
    os.environ['OMP_NUM_THREADS'] = '4'  # Limit parallelism
    os.environ['MKL_NUM_THREADS'] = '1'

    # For compatibility
    os.environ['GFORTRAN_UNBUFFERED_ALL'] = '1'
```

Fallback for large molecules:
```python
def generate_conformers_adaptive(mol: Chem.Mol):
    """Choose conformer method based on molecule size."""
    n_atoms = mol.GetNumAtoms()

    if n_atoms > 500:
        logger.warning(f"Molecule has {n_atoms} atoms, skipping CREST (xTB memory issues)")
        return generate_conformers_rdkit(mol, num_confs=30)
    elif n_atoms > 200:
        logger.info(f"Large molecule ({n_atoms} atoms), using RDKit instead of CREST")
        return generate_conformers_rdkit(mol, num_confs=50)
    else:
        return generate_conformers_crest(mol)
```

**Warning signs:**
- Process killed with signal 9 (SIGKILL)
- No error message in CREST/xTB output (killed externally)
- `dmesg` shows "Out of memory: Killed process"
- CREST stops mid-run on large molecules

---

### Gotcha 2: CREST ALPB vs GBSA Solvation Model Confusion

**Problem:**
Passing `--gbsa chcl3` to CREST fails silently or gives poor results. GBSA is parameterized for GFN1-xTB and GFN2-xTB, but not all solvents. ALPB is newer, more accurate, but has different solvent name syntax.

**Root cause:**
- ALPB is "more recent development offering improved accuracy over GBSA" ([xTB Docs](https://xtb-docs.readthedocs.io/en/latest/gbsa.html))
- GBSA estimates electrostatic contribution using generalized Born, ALPB solves linearized Poisson-Boltzmann analytically
- Different solvent compatibility: "ALPB model is parameterized for GFN1-xTB, GFN2-xTB, and GFN-FF, but not for GFN0-xTB" ([Grimme-Lab Workshop](https://grimme-lab.github.io/workshops/page/xtb/Solvation))
- Solvent names differ: ALPB uses `water`, GBSA uses full solvent database

**Prevention:**
Document solvation model choice explicitly:
```python
CREST_SOLVENTS = {
    'water': 'water',
    'chcl3': 'chcl3',
    'dmso': 'dmso',
    'acetone': 'acetone',
    'methanol': 'methanol',
}

def generate_conformers_crest(mol, solvent: str = None):
    """Generate conformers with CREST using ALPB solvation."""
    cmd = ['crest', xyz_file, '--gfn2']

    if solvent and solvent in CREST_SOLVENTS:
        # Use ALPB (more accurate than GBSA)
        cmd.extend(['--alpb', CREST_SOLVENTS[solvent]])
    elif solvent:
        raise ValueError(
            f"Solvent '{solvent}' not supported by CREST ALPB. "
            f"Supported: {list(CREST_SOLVENTS.keys())}"
        )

    # Run CREST...
```

**Warning signs:**
- CREST outputs "WARNING: Unknown solvent for ALPB"
- Conformer energies differ drastically between CREST (gas) and NWChem (COSMO)
- Hydration energies unrealistic (>10 kcal/mol error vs experimental)

---

### Gotcha 3: RDKit MMFF Optimization Failures on Macrocycles

**Problem:**
`MMFFOptimizeMoleculeConfs()` returns `not_converged=1` for macrocycles (rings >12 atoms), strained rings, or molecules with complex stereochemistry. Conformers have unrealistic geometries (stretched bonds, non-planar aromatic rings). Referenced research notes "macrocycles... create a challenge for conformer generators" due to "larger number of degrees of freedom" ([RDKit Blog](https://greglandrum.github.io/rdkit-blog/posts/2023-05-17-understanding-confgen-errors.html)).

**Root cause:**
- MMFF94 force field struggles with ring strain
- Default `maxIters=200` insufficient for large flexible systems
- `EmbedMultipleConfs` can give "conformations that are not physically realistic (stretched saturated rings, non-linear sp hybridized carbons)" ([RDKit Discussion #5368](https://github.com/rdkit/rdkit/discussions/5368))

**Prevention:**
Check convergence and filter failed conformers:
```python
def generate_conformers_rdkit_robust(mol: Chem.Mol, num_confs: int = 50):
    """Generate conformers with convergence checking."""
    mol = Chem.AddHs(mol)

    params = AllChem.KDGParams()
    params.pruneRmsThresh = 0.5
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    # MMFF optimization with increased max iterations for macrocycles
    results = AllChem.MMFFOptimizeMoleculeConfs(
        mol,
        numThreads=0,
        maxIters=500,  # Increase from default 200
    )

    # Filter to converged conformers only
    converged_cids = []
    for i, (not_converged, energy) in enumerate(results):
        if not_converged == 0:  # 0 = converged
            converged_cids.append(cids[i])
        else:
            logger.warning(f"Conformer {i} MMFF optimization did not converge")

    if len(converged_cids) == 0:
        raise RuntimeError(
            "All MMFF optimizations failed. Consider using CREST instead."
        )

    if len(converged_cids) < num_confs * 0.5:
        logger.warning(
            f"Only {len(converged_cids)}/{num_confs} conformers converged. "
            "Molecule may be too complex for MMFF."
        )

    return mol, converged_cids
```

**Warning signs:**
- `not_converged=1` for >50% of conformers
- Conformer energies extremely high (>100 kcal/mol)
- Visual inspection shows distorted geometries
- DFT optimization drastically changes geometry (RMSD >2 Å vs MMFF)

---

### Gotcha 4: RDKit Stereochemistry Loss During Conformer Generation

**Problem:**
Generating conformers from SMILES with defined stereochemistry (`CC[C@H](O)C`), but output conformers have wrong or unspecified stereochemistry. Averaging over enantiomers instead of single stereoisomer.

**Root cause:**
- Issue documented: "Stereochemistry doesn't retain in the output structures after the generation of conformers" ([GitHub Issue #1504](https://github.com/rdkit/rdkit/issues/1504))
- Conformer generation checks fail with `BAD_DOUBLE_BOND_STEREO` when stereo constraints can't be satisfied ([RDKit Blog](https://greglandrum.github.io/rdkit-blog/posts/2023-05-17-understanding-confgen-errors.html))
- Distance geometry can generate geometries that violate stereochemistry

**Prevention:**
Preserve stereochemistry explicitly:
```python
def generate_conformers_preserve_stereo(mol: Chem.Mol, num_confs: int = 50):
    """Generate conformers while preserving stereochemistry."""
    # Ensure stereochemistry is defined
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    mol = Chem.AddHs(mol)

    params = AllChem.KDGParams()
    params.pruneRmsThresh = 0.5
    params.useRandomCoords = False  # Use 3D coords if available

    # Track failures
    params.trackFailures = True  # Available in RDKit 2023.03.1+

    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    if len(cids) == 0:
        # Check what failed
        failures = params.getFailures()  # Dict of failure types
        if 'BAD_DOUBLE_BOND_STEREO' in failures:
            raise ValueError(
                "Conformer generation failed due to stereochemistry constraints. "
                "Check SMILES stereochemistry specification."
            )

    # Verify stereochemistry preserved
    for cid in cids:
        Chem.AssignStereochemistry(mol, confId=cid, force=True)
        # Could check against original stereochemistry here

    return mol, cids
```

**Warning signs:**
- NMR predictions symmetric when experimental spectra are asymmetric
- Generated conformers include both R and S enantiomers
- RDKit logs `BAD_DOUBLE_BOND_STEREO` failures
- Conformer count much lower than requested

---

## Performance Traps

### Trap 1: Sequential Conformer DFT Optimization

**Problem:**
Running 50 conformer DFT optimizations sequentially takes 50 × 15min = 12.5 hours. Single Huey worker processes one conformer at a time, wasting CPU during I/O-bound NWChem phases.

**Why it happens:**
- Current architecture: single Huey worker on single VM ([PROJECT.md](https://github.com))
- Huey task runs entire conformer ensemble as one job
- NWChem launched sequentially in Python loop

**Mitigation (within single-worker constraint):**
Use subprocess parallelism for independent conformer optimizations:
```python
import concurrent.futures
from pathlib import Path

def optimize_conformers_parallel(
    conformers: list[ConformerData],
    job_dir: Path,
    preset: dict,
    solvent: str,
    max_workers: int = 4,  # CPU cores / processes_per_nwchem
):
    """Optimize conformers in parallel using subprocess pool."""

    def optimize_one(conf_id: int, conf_data: ConformerData):
        return run_calculation(
            smiles=conf_data.smiles,
            job_dir=job_dir,
            conformer_id=conf_id,  # Unique scratch dir
            preset=preset,
            solvent=solvent,
            processes=1,  # Low MPI parallelism per job
        )

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(optimize_one, i, conf): i
            for i, conf in enumerate(conformers)
        }

        results = {}
        for future in concurrent.futures.as_completed(futures):
            conf_id = futures[future]
            try:
                results[conf_id] = future.result()
            except Exception as e:
                logger.error(f"Conformer {conf_id} optimization failed: {e}")
                # Continue with other conformers

        return results
```

**Trade-off:**
- Parallel: 50 conformers / 4 workers = 12.5 rounds × 15min = ~3 hours (4× speedup)
- Must reduce NWChem MPI processes from 4 to 1 per job (or use `max_workers=1`)
- Scratch directory isolation (Pitfall 5) becomes critical

**Warning signs:**
- CPU usage <25% during conformer optimizations (underutilization)
- Wall time scales linearly with conformer count (no parallelism)
- Single NWChem process uses <100% CPU (serial execution)

---

### Trap 2: Disk Space Explosion from Conformer NWChem Calculations

**Problem:**
50 conformers × 2 NWChem jobs (opt + NMR) × 500 MB scratch = 50 GB disk usage. Scratch files not cleaned up between conformers. Disk fills, job fails mid-run.

**Why it happens:**
- NWChem generates large scratch files (`.db`, `.movecs`, `.gridpts`)
- Current code does not clean scratch between steps
- Conformer ensemble multiplies disk usage by N

**Prevention:**
Aggressive cleanup after each conformer:
```python
def run_calculation_with_cleanup(conformer_id: int, ...):
    """Run NWChem calculation with immediate cleanup."""
    scratch_dir = job_dir / "scratch" / f"conf_{conformer_id}"
    scratch_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Run optimization
        result_opt = run_nwchem(opt_input, opt_output, processes=4)

        # Extract geometry immediately
        geometry = extract_optimized_geometry(opt_output.read_text())

        # Clean optimization scratch files (keep only .out for debugging)
        for pattern in ['*.db', '*.movecs', '*.gridpts', '*.zmat']:
            for f in scratch_dir.glob(pattern):
                f.unlink()

        # Run NMR calculation
        result_nmr = run_nwchem(nmr_input, nmr_output, processes=4)

        # Extract NMR data immediately
        nmr_data = parse_shielding_output(nmr_output.read_text())

        # Clean NMR scratch files
        for pattern in ['*.db', '*.movecs', '*.gridpts', '*.zmat']:
            for f in scratch_dir.glob(pattern):
                f.unlink()

        return {'geometry': geometry, 'nmr': nmr_data}

    finally:
        # Final cleanup (keep .out files for debugging)
        pass
```

Monitor disk usage:
```python
import shutil

def check_disk_space(job_dir: Path, required_gb: float = 10.0):
    """Check available disk space before starting calculation."""
    stat = shutil.disk_usage(job_dir)
    available_gb = stat.free / (1024**3)

    if available_gb < required_gb:
        raise RuntimeError(
            f"Insufficient disk space: {available_gb:.1f} GB available, "
            f"{required_gb:.1f} GB required"
        )
```

**Warning signs:**
- Job fails with "No space left on device"
- `du -sh data/jobs/{job_id}` shows >10 GB per job
- Scratch directory contains hundreds of `.db` files
- Disk usage increases monotonically during conformer loop

---

### Trap 3: Long-Running Ensemble Jobs Timing Out

**Problem:**
Web service has 30-minute request timeout. Conformer ensemble job takes 3 hours. Job starts, timeout fires, user sees error, but job continues running in background. User re-submits, creating duplicate jobs.

**Why it happens:**
- HTTP request timeout != background job timeout
- Huey tasks run asynchronously, independent of web request
- Current status.json doesn't indicate "job is still running" after request times out

**Prevention:**
Educate users with realistic time estimates:
```python
def estimate_job_time(smiles: str, use_crest: bool, num_conformers: int = None) -> int:
    """Estimate job completion time in minutes."""
    mol = Chem.MolFromSmiles(smiles)
    n_atoms = mol.GetNumAtoms()
    n_rotatable = Chem.Descriptors.NumRotatableBonds(mol)

    if not use_crest and num_conformers is None:
        num_conformers = calculate_conformer_count(mol)

    # Single-conformer baseline: ~5 min opt + ~10 min NMR
    single_conf_time = 5 + (n_atoms / 10) + 10 + (n_atoms / 5)

    if num_conformers == 1:
        return int(single_conf_time)

    # Conformer generation time
    if use_crest:
        crest_time = 30 + (n_rotatable * 5)  # 30-90 min
    else:
        crest_time = 1  # RDKit fast

    # DFT optimization: parallel speedup if max_workers > 1
    opt_time = num_conformers * single_conf_time / 4  # 4 workers

    total_minutes = crest_time + opt_time
    return int(total_minutes * 1.2)  # Add 20% buffer
```

Show estimate in UI before submission:
```python
# In API endpoint
estimated_minutes = estimate_job_time(smiles, use_crest=True)
return {
    "job_id": job_id,
    "status": "queued",
    "estimated_completion_minutes": estimated_minutes,
    "message": f"Job queued. Estimated time: {estimated_minutes} minutes. "
               f"You will receive an email when complete."
}
```

**Warning signs:**
- Multiple jobs with identical SMILES from same user
- Jobs stuck in "running" state for >2 hours
- Users reporting "timeout" but job completes later
- Status endpoint returns 504 Gateway Timeout

---

## "Looks Done But Isn't" Checklist

### Missing Feature 1: Partial Failure Handling

**Scenario:**
Ensemble has 50 conformers. 45 DFT optimizations converge, 5 fail (SCF convergence issues). Current code crashes entire job on first failure, discarding 45 successful results.

**What it should do:**
Continue with successful conformers, log failures, adjust Boltzmann weights to sum of successful conformers only.

```python
def optimize_conformers_robust(conformers: list, ...):
    """Optimize conformers with partial failure handling."""
    successful = []
    failed = []

    for i, conf in enumerate(conformers):
        try:
            result = run_calculation(conf, ...)
            successful.append({'id': i, 'result': result})
        except RuntimeError as e:
            logger.warning(f"Conformer {i} optimization failed: {e}")
            failed.append({'id': i, 'error': str(e)})

    if len(successful) == 0:
        raise RuntimeError("All conformer optimizations failed")

    if len(failed) > 0:
        logger.warning(
            f"{len(failed)}/{len(conformers)} conformers failed. "
            f"Proceeding with {len(successful)} successful conformers."
        )

    return successful, failed
```

**Detection:**
- Log shows "X conformers failed" but job status is "complete"
- nmr_results.json includes metadata: `"conformers_attempted": 50, "conformers_successful": 45`

---

### Missing Feature 2: Progress Reporting for Multi-Step Pipeline

**Scenario:**
User sees "Running" for 2 hours with no updates. Are conformers being generated? Is DFT running? Which conformer (5/50)? No visibility into pipeline progress.

**What it should do:**
Update `current_step` and `current_step_label` for each major stage:

```python
# In conformer ensemble task
start_step(job_id, "conformer_generation", "Generating conformers with CREST")
conformers = generate_conformers_crest(...)

start_step(job_id, "conformer_filtering", "Filtering conformers by energy")
conformers = filter_by_energy(conformers, window=6.0)

for i, conf in enumerate(conformers):
    start_step(
        job_id,
        f"conformer_opt_{i}",
        f"Optimizing conformer {i+1}/{len(conformers)}"
    )
    optimize_conformer(conf)

start_step(job_id, "conformer_nmr", "Computing NMR shieldings for ensemble")
# ... NMR calculations

start_step(job_id, "boltzmann_averaging", "Averaging NMR shifts")
# ... averaging
```

**Detection:**
- Frontend can display: "Optimizing conformer 23/50 (Estimated: 45 min remaining)"
- Steps are logged in `steps_completed` for performance analysis

---

### Missing Feature 3: Degeneracy Handling for Single Conformer

**Scenario:**
Rigid molecule (e.g., benzene) has only 1 conformer. Code calculates `weights = [1.0]` and averages over single value, producing correct result. But Boltzmann formula `exp(-0/RT) / sum(exp(-0/RT)) = 1/1 = 1` is degenerate case.

**What it should do:**
Detect single-conformer case and skip averaging:

```python
def boltzmann_average_shifts(conformers: list, energies: np.ndarray):
    """Average NMR shifts with Boltzmann weighting."""
    if len(conformers) == 1:
        logger.info("Single conformer, no averaging needed")
        return conformers[0]['nmr_shifts']

    weights = boltzmann_weights(energies)
    # ... averaging
```

**Why it matters:**
- Clarity: Log shows "no averaging needed" vs silent averaging
- Performance: Skip unnecessary computation
- Numerical stability: Avoid `1/1` division

---

### Missing Feature 4: Symmetry-Equivalent Conformer Overcounting

**Scenario:**
n-butane has 2 gauche conformers that are mirror images (g+ and g-). CREST generates both, each with weight 0.25. Averaging counts them separately, but they should have combined degeneracy factor of 2 (physical reality: ggauche = 2, gtrans = 1).

**What it should do:**
Detect symmetry-equivalent conformers and apply degeneracy factors:

```python
def calculate_degeneracy_factors(conformers: list) -> list[int]:
    """Calculate symmetry degeneracy for each conformer.

    Conformers that are mirror images (enantiomers) should have degeneracy = 2
    if both are present, or degeneracy = 1 if only one is kept.
    """
    # This is complex - requires symmetry detection
    # Simplified: assume CREST deduplication handles this
    # Advanced: use RDKit symmetry perception or point group analysis

    return [1] * len(conformers)  # Default: no correction
```

Reference notes: "Conformers can come in different absolute configurations (mirror images), and this affects the degeneracy factor (g)" ([Exploring Conformational Averaging](https://pmc.ncbi.nlm.nih.gov/articles/PMC8139173/)).

**Detection:**
- Pairwise RMSD shows near-identical conformers (RMSD < 0.1 Å after reflection)
- Manual inspection: two conformers differ only by chirality flip
- Literature comparison: experimental NMR matches better when degeneracy applied

---

## Pitfall-to-Phase Mapping

| Phase | Critical Pitfalls to Address | Technical Debt to Address |
|-------|------------------------------|---------------------------|
| **Phase 1: Conformer Data Model** | Pitfall 4 (Atom ordering) | Pattern 1 (Adaptive conformer count) |
| **Phase 2: RDKit Conformer Generation** | None (RDKit is simpler) | Pattern 2 (Deduplication), Pattern 3 (KDG vs ETKDG) |
| **Phase 3: DFT Multi-Conformer Optimization** | Pitfall 5 (Scratch directory conflicts), Pitfall 6 (Energy source) | Trap 1 (Sequential optimization), Trap 2 (Disk space) |
| **Phase 4: CREST Integration** | Pitfall 3 (CREST hanging), Gotcha 1 (xTB memory), Gotcha 2 (ALPB vs GBSA) | None |
| **Phase 5: Boltzmann Averaging** | Pitfall 1 (Numerical overflow), Pitfall 2 (Energy units), Pitfall 4 (Atom ordering) | Missing Feature 1 (Partial failures), Missing Feature 4 (Symmetry) |
| **Phase 6: Integration & UI** | None | Trap 3 (Timeout handling), Missing Feature 2 (Progress reporting) |

---

## Sources

### CREST/xTB Issues
- [CREST Issue #242: Cregen eliminates valid conformers](https://github.com/crest-lab/crest/issues/242)
- [CREST Issue #105: Initial geometry optimization failed](https://github.com/crest-lab/crest/issues/105)
- [CREST Issue #313: Initial geometry optimization failed](https://github.com/crest-lab/crest/issues/313)
- [CREST Issue #382: Conformational sampling help](https://github.com/crest-lab/crest/issues/382)
- [xTB Issue #700: Crash with large systems](https://github.com/grimme-lab/xtb/issues/700)
- [xTB Discussion #811: xTB killed with larger molecules](https://github.com/grimme-lab/xtb/discussions/811)
- [xTB Issue #439: XTB crashes with large molecules](https://github.com/grimme-lab/xtb/issues/439)
- [xTB Issue #242: gfnff crashes with large molecules](https://github.com/grimme-lab/xtb/issues/242)

### RDKit Conformer Generation
- [RDKit Issue #5250: Conformer Generation takes very long time](https://github.com/rdkit/rdkit/issues/5250)
- [RDKit Blog: Understanding conformer generation failures](https://greglandrum.github.io/rdkit-blog/posts/2023-05-17-understanding-confgen-errors.html)
- [RDKit Discussion #8226: Best practices - Conformer generation](https://github.com/rdkit/rdkit/discussions/8226)
- [RDKit Issue #8001: RMS pruning misses conformers](https://github.com/rdkit/rdkit/issues/8001)
- [RDKit Issue #1504: Stereochemistry doesn't retain](https://github.com/rdkit/rdkit/issues/1504)
- [RDKit Discussion #5368: Constrained optimization](https://github.com/rdkit/rdkit/discussions/5368)

### Boltzmann Weighting
- [Evaluating Error in Boltzmann Weighting](https://corinwagen.github.io/public/blog/20221228_boltzmann_error.html)
- [Exp-normalize trick — Graduate Descent](https://timvieira.github.io/blog/post/2014/02/11/exp-normalize-trick/)
- [Exploring the impacts of conformer selection methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8139173/)

### ALPB/GBSA Solvation
- [xTB Docs: Implicit Solvation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html)
- [Robust and Efficient Implicit Solvation Model](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.1c00471)
- [Grimme-Lab Workshop: Solvation](https://grimme-lab.github.io/workshops/page/xtb/Solvation)

### Atom Ordering and Stereochemistry
- [Nature Communications 2025: DFT calculations and theory do not support enantiospecificity](https://www.nature.com/articles/s41467-025-66247-0)
- [May the Force (Field) Be with You: Conformational Searches](https://pmc.ncbi.nlm.nih.gov/articles/PMC9694776/)

### NWChem Scratch
- [NWChem Wiki: Scratch_Dir](https://github.com/nwchemgit/nwchem/wiki/Scratch_Dir)

### Energy Units
- [Energy Units Converter - Hacettepe](https://yunus.hacettepe.edu.tr/~ttugsuz/Hartree.html)
- [Energy Conversion Table](https://wild.life.nctu.edu.tw/class/common/energy-unit-conv-table.html)

---

*Pitfalls research for: Conformational Sampling NMR*
*Researched: 2026-01-26*
*Confidence: HIGH (based on GitHub issues, official documentation, and computational chemistry literature)*
