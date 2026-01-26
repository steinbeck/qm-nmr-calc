# Stack Research: Conformational Sampling (v2.0)

**Domain:** Conformational Sampling with Boltzmann-weighted NMR Averaging
**Milestone:** v2.0 - Add conformational sampling features
**Researched:** 2026-01-26
**Confidence:** HIGH

## Executive Summary

Adding conformational sampling to qm-nmr-calc requires **minimal new dependencies**. The core computational libraries (scipy with numpy, RDKit) are already present. The primary additions are optional external binaries (CREST/xTB) for high-accuracy conformer generation, which should be detected on PATH but not required for basic functionality.

**Key insight:** RDKit's KDG (Knowledge Distance Geometry) method is ideal for solution-phase NMR conformer generation because it avoids the crystal structure bias inherent in ETKDG.

**Architecture decision:** Two-tier approach with RDKit KDG as the always-available baseline and CREST as an optional high-accuracy upgrade.

## Stack Status: Existing vs New

### Already Available (No Changes Needed)

| Technology | Current Version | Purpose for v2.0 |
|------------|----------------|------------------|
| **RDKit** | >=2025.09.3 | Conformer generation using KDG method |
| **scipy** | >=1.17.0 | Boltzmann statistics (includes numpy) |
| **numpy** | (via scipy) | Array operations, exponential functions for weighting |
| **NWChem** | >=7.2.0 | DFT optimization and NMR for each conformer |
| **subprocess** | stdlib | Calling CREST/xTB binaries |
| **pathlib** | stdlib | File path handling for conformer ensembles |

### New Optional Dependencies (v2.0)

| Technology | Version | Purpose | Required? |
|------------|---------|---------|-----------|
| **CREST binary** | 3.0.2 or continuous | Metadynamics-based conformer generation | Optional (auto-detected on PATH) |
| **xTB binary** | 6.7.1 | GFN2-xTB (required by CREST for some features) | Optional (with CREST) |

### What We're NOT Adding

- ❌ **xtb-python** - Not needed, call CREST binary directly
- ❌ **OpenBabel** - Already an ISiCLE dependency for format conversion
- ❌ **CENSO** - Excellent tool but adds ORCA/Turbomole dependency we don't need
- ❌ **scipy.stats.boltzmann** - Wrong distribution (PMF), use custom numpy implementation
- ❌ **New Python packages** - Zero new Python dependencies needed

## Detailed Stack Components

### 1. RDKit KDG (Core Conformer Generation)

**Status:** Already available (rdkit>=2025.09.3 in pyproject.toml)

**What is KDG?**
- **K**nowledge Distance **G**eometry - uses basic chemical knowledge without crystal bias
- Developed by Riniker & Landrum (2015), available in RDKit since 2015.09.1
- Ideal for solution-phase conformers (unlike ETKDG which has crystal structure bias from CSD)

**Why KDG over ETKDG for NMR?**

From RDKit community discussion:
> "The ET terms bias the sampling of conformer space to favor torsion values that have been observed in crystal structures. These are quite useful if you want conformers for condensed phases or for protein-ligand docking, but can miss regions of conformer space that you may want to sample if you are working with compounds in the gas phase or solution. For those cases you are probably best served by sticking with KDG, not ETKDGv3."

**API Usage:**
```python
from rdkit.Chem import AllChem

# KDG for unbiased solution-phase sampling
params = AllChem.KDG()  # Returns EmbedParameters object
params.pruneRmsThresh = 0.5  # RMSD-based conformer filtering
params.numThreads = 0  # Use all available CPU cores

# Generate conformers
cids = AllChem.EmbedMultipleConfs(mol, numConfs=50, params=params)

# MMFF optimization
results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
```

**Version Timeline:**
- 2015.09.1: KDG introduced
- 2024.03: ETKDGv3 became default (not suitable for solution NMR)
- 2025.09.1: Allene/cumulene support added to KDG
- 2025.09.3: Current project dependency (confirmed in pyproject.toml)

**Confidence:** HIGH - Verified in official RDKit documentation, confirmed available in project dependencies

### 2. scipy/numpy (Boltzmann Weighting)

**Status:** Already available (scipy>=1.17.0 in pyproject.toml)

**What we use it for:**
- `numpy.exp()` for Boltzmann factor calculations: `exp(-ΔE/RT)`
- `numpy.sum()` for normalization
- `numpy.min()` for energy referencing
- Array operations for vectorized calculations

**Why NOT scipy.stats.boltzmann?**

scipy.stats.boltzmann is for discrete Boltzmann distributions (probability mass function), not conformer weighting. We need a custom implementation:

```python
import numpy as np

def boltzmann_weights(energies_kcal, temperature=298.15):
    """Calculate Boltzmann weights from relative energies."""
    R = 0.001987204  # kcal/(mol·K)
    RT = R * temperature

    # Reference to minimum energy
    rel_energies = energies_kcal - np.min(energies_kcal)

    # Boltzmann factors
    boltz_factors = np.exp(-rel_energies / RT)

    # Normalize to sum = 1.0
    weights = boltz_factors / np.sum(boltz_factors)

    return weights
```

**Confidence:** HIGH - scipy confirmed in pyproject.toml, standard statistical mechanics implementation

### 3. CREST Binary (Optional High-Accuracy Mode)

**Status:** NEW optional dependency (auto-detected on system PATH)

**Latest Version:**
- Stable: 3.0.2 (August 25, 2020)
- Continuous: December 23, 2025 (pre-release, minimal testing)
- **Recommended:** Install from conda-forge (automatically gets latest stable)

**What is CREST?**
- **C**onformer-**R**otamer **E**nsemble **S**ampling **T**ool
- Developed by Grimme lab (University of Bonn)
- Uses metadynamics with GFN2-xTB for conformer exploration
- Published 2020, major update 2024 with refactored backend

**Why CREST?**
- State-of-the-art conformer sampling via metadynamics
- GFN2-xTB energies have better correlation with DFT than force fields (Spearman ρ ~ 0.39-0.47 vs -0.1 to -0.45 for MMFF)
- Handles complex molecules (macrocycles, highly flexible systems)
- Solvation models available (ALPB, GBSA)

**Installation:**

```bash
# Recommended: conda-forge (easiest)
conda install -c conda-forge crest xtb

# Alternative: Download binary from GitHub
# https://github.com/crest-lab/crest/releases
```

**Command Line Interface:**

```bash
crest input.xyz \
    --gfn2 \                    # Use GFN2-xTB (default, best)
    --chrg 0 \                  # Molecular charge
    --ewin 6.0 \                # Energy window in kcal/mol (default 6.0)
    --alpb chcl3 \              # ALPB implicit solvation (better than GBSA)
    --temp 298.15 \             # Temperature for Boltzmann populations
    --T 8                       # Number of CPU threads
```

**Output Format:**

CREST writes `crest_conformers.xyz` in multi-structure XYZ format:

```
[num_atoms]
[energy_in_Hartree]  [optional metadata]
Element1  x1  y1  z1
Element2  x2  y2  z2
...
[num_atoms]
[energy_in_Hartree]
Element1  x1  y1  z1
...
```

**Key CLI Options:**

| Option | Default | Purpose |
|--------|---------|---------|
| `--gfn2` | yes | Use GFN2-xTB method |
| `--ewin <real>` | 6.0 | Energy window (kcal/mol) |
| `--alpb <solvent>` | none | ALPB implicit solvation |
| `--gbsa <solvent>` | none | GBSA implicit solvation |
| `--chrg <int>` | 0 | Molecular charge |
| `--temp <real>` | 298.15 | Temperature (K) for Boltzmann |
| `--T <int>` | auto | Number of CPU threads |

**Integration Pattern:**

```python
import subprocess
import shutil

def has_crest():
    """Check if CREST is available on PATH."""
    return shutil.which('crest') is not None

def generate_conformers(mol, method='auto'):
    """Generate conformers using best available method."""
    if method == 'crest' and has_crest():
        return generate_conformers_crest(mol)
    elif method == 'auto' and has_crest():
        # Use CREST for flexible molecules (>5 rotatable bonds)
        n_rotatable = Descriptors.NumRotatableBonds(mol)
        if n_rotatable > 5:
            return generate_conformers_crest(mol)

    # Fallback to RDKit KDG (always available)
    return generate_conformers_rdkit_kdg(mol)
```

**Confidence:** HIGH - Versions verified from GitHub releases, CLI options from official documentation

### 4. xTB Binary (Required by CREST)

**Status:** NEW optional dependency (installed with CREST)

**Latest Version:** 6.7.1 (July 23, 2024)

**Relationship to CREST:**

From CREST documentation:
> "Most CREST applications will require access to the xtb program, and you will need to install xtb on your machine and make sure it can be executed globally."

However, for CREST 3.0+:
> "While xtb is technically not needed for the primary runtypes of CREST versions >3.0 thanks to an integration of tblite, some functionalities, like QCG, still require it!"

**Recommendation:** Always install xTB alongside CREST. Conda handles this automatically:

```bash
conda install -c conda-forge crest xtb
```

**What is GFN2-xTB?**
- **G**eometry, **F**requency, and **N**on-covalent interactions **xTB** method
- Semi-empirical quantum method (faster than DFT)
- Parametrized for elements up to Z=86
- Accuracy: ~2 kcal/mol MAE for conformer relative energies

**Why GFN2-xTB energies are good for ranking:**

Research shows GFN2-xTB conformer energies correlate reasonably with DFT (Spearman ρ ~ 0.39-0.47), whereas MMFF can anticorrelate (ρ ~ -0.1 to -0.45). This makes xTB suitable for pre-DFT filtering but not final Boltzmann weighting.

**Confidence:** HIGH - Version verified from GitHub releases, relationship to CREST confirmed in official docs

## CREST Output Parsing

**No external library needed.** Simple Python parser using stdlib:

```python
def parse_crest_ensemble(filepath):
    """
    Parse CREST multi-xyz output file.

    Args:
        filepath: Path to crest_conformers.xyz

    Returns:
        List of dicts with 'n_atoms', 'energy_hartree', 'xyz'
    """
    conformers = []

    with open(filepath) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        # Read atom count
        n_atoms = int(lines[i].strip())
        i += 1

        # Read energy from comment line (first token in Hartree)
        energy_line = lines[i].strip().split()
        energy_hartree = float(energy_line[0])
        i += 1

        # Read coordinate lines
        coords_lines = []
        for _ in range(n_atoms):
            coords_lines.append(lines[i].strip())
            i += 1

        conformers.append({
            'n_atoms': n_atoms,
            'energy_hartree': energy_hartree,
            'xyz': '\n'.join(coords_lines)
        })

    return conformers
```

**Format notes:**
- Structures must have same atom count
- Atom order must not change between conformers
- Energy in Hartree (1 Hartree = 627.509 kcal/mol)
- Structures pre-sorted by energy (lowest first)

**Confidence:** HIGH - Format documented in official CREST documentation

## Boltzmann Weighting Implementation

**No external library needed.** Simple numpy implementation:

```python
import numpy as np

def boltzmann_weights(energies_kcal, temperature=298.15):
    """
    Calculate Boltzmann weights from relative energies.

    Args:
        energies_kcal: Array of conformer energies in kcal/mol
        temperature: Temperature in Kelvin (default 298.15 K)

    Returns:
        Array of normalized Boltzmann weights (sum = 1.0)

    Notes:
        - Input energies can be absolute or relative
        - Function references to minimum energy automatically
        - Uses gas constant R = 0.001987204 kcal/(mol·K)
    """
    R = 0.001987204  # kcal/(mol·K) - gas constant
    RT = R * temperature

    # Reference to minimum energy (ensure ΔE_min = 0)
    rel_energies = energies_kcal - np.min(energies_kcal)

    # Calculate Boltzmann factors: exp(-ΔE/RT)
    boltz_factors = np.exp(-rel_energies / RT)

    # Normalize to sum = 1.0
    weights = boltz_factors / np.sum(boltz_factors)

    return weights

def boltzmann_average(values, weights):
    """
    Calculate Boltzmann-weighted average.

    Args:
        values: Array of property values per conformer
        weights: Array of Boltzmann weights (must sum to 1.0)

    Returns:
        Weighted average: Σ(w_i × v_i)
    """
    return np.sum(weights * values)

def filter_by_population(energies_kcal, threshold=0.01, temperature=298.15):
    """
    Filter conformers by Boltzmann population threshold.

    Args:
        energies_kcal: Array of conformer energies
        threshold: Minimum population (default 1%)
        temperature: Temperature in K

    Returns:
        Boolean array indicating conformers to keep
    """
    weights = boltzmann_weights(energies_kcal, temperature)
    return weights >= threshold

def filter_by_cumulative(energies_kcal, cumulative=0.95, temperature=298.15):
    """
    Filter conformers to capture target cumulative population.

    Args:
        energies_kcal: Array of conformer energies
        cumulative: Target cumulative population (default 95%)
        temperature: Temperature in K

    Returns:
        Boolean array indicating conformers to keep

    Notes:
        - Conformers must be sorted by energy (lowest first)
        - Keeps conformers until cumulative population reached
    """
    weights = boltzmann_weights(energies_kcal, temperature)
    cumsum = np.cumsum(weights)
    return cumsum <= cumulative
```

**Confidence:** HIGH - Standard statistical mechanics, implementation verified against reference document

## Energy Windows and Filtering Strategy

### Pre-DFT Filtering (After CREST/RDKit)

**Energy window:** 5-6 kcal/mol above global minimum

**Rationale:**
- At 298 K, conformer at 6 kcal/mol contributes <0.01% to population
- Captures >99% of thermal population
- Avoids wasting DFT resources on negligible contributors

**Implementation:**
```python
def filter_pre_dft(conformers, energy_window=6.0):
    """Filter conformers before expensive DFT calculations."""
    energies = np.array([c['energy_kcal'] for c in conformers])
    min_energy = np.min(energies)
    keep = (energies - min_energy) <= energy_window
    return [c for c, k in zip(conformers, keep) if k]
```

### Post-DFT Filtering (After Geometry Optimization)

**Strategy:** Keep conformers representing top 95-99% cumulative Boltzmann population

**Rationale:**
- DFT energies are more accurate than xTB/MMFF
- Energy ordering may change after DFT optimization
- Adaptive threshold ensures meaningful contributors retained
- 3 kcal/mol window captures similar population but may miss edge cases

**Implementation:**
```python
def filter_post_dft(conformers, cumulative=0.95, temperature=298.15):
    """Filter conformers after DFT optimization."""
    # Sort by DFT energy
    conformers = sorted(conformers, key=lambda c: c['dft_energy_kcal'])
    energies = np.array([c['dft_energy_kcal'] for c in conformers])

    # Calculate cumulative population
    weights = boltzmann_weights(energies, temperature)
    cumsum = np.cumsum(weights)

    # Keep conformers until target cumulative reached
    keep_count = np.sum(cumsum <= cumulative) + 1  # +1 to include threshold-crossing conformer

    return conformers[:keep_count]
```

### Summary Table

| Stage | Filter Type | Threshold | Rationale |
|-------|------------|-----------|-----------|
| **Pre-DFT** | Energy window | 6 kcal/mol | Captures >99% population, reduces DFT cost |
| **Post-DFT** | Cumulative population | 95-99% | Adaptive, accounts for DFT reranking |
| **Alternative** | Population threshold | >1% | Fixed threshold, simpler but less robust |

**Confidence:** HIGH - Based on thermal population statistics and verified against reference document

## Conformer Count Guidelines

| Molecule Type | Rotatable Bonds | Method | Initial Conformers | Expected After Filtering |
|---------------|-----------------|--------|-------------------|-------------------------|
| **Rigid** | 0-1 | RDKit KDG | 5-10 | 1-3 |
| **Semi-rigid** | 2-4 | RDKit KDG | 20-50 | 5-15 |
| **Flexible** | 5-8 | CREST or RDKit | 50-100 | 10-30 |
| **Very flexible** | >8 | CREST required | 100-500 | 20-50 |

**Adaptive strategy:**
```python
from rdkit.Chem import Descriptors

def recommend_conformer_count(mol):
    """Recommend initial conformer count based on flexibility."""
    n_rotatable = Descriptors.NumRotatableBonds(mol)

    if n_rotatable <= 1:
        return 10
    elif n_rotatable <= 4:
        return 30
    elif n_rotatable <= 8:
        return 100
    else:
        return 200

def recommend_method(mol):
    """Recommend conformer generation method."""
    n_rotatable = Descriptors.NumRotatableBonds(mol)

    if n_rotatable <= 4:
        return 'rdkit'  # KDG sufficient
    elif n_rotatable <= 8 and has_crest():
        return 'crest'  # CREST preferred if available
    else:
        return 'rdkit'  # Fallback (extensive sampling)
```

**Confidence:** HIGH - Guidelines from conformer generation literature and reference document

## Computational Cost Estimation

### Per-conformer costs (typical 30-atom drug-like molecule)

| Stage | Time | Parallelizable | Notes |
|-------|------|----------------|-------|
| RDKit KDG generation | ~1 sec | Yes (numThreads) | 50 conformers in ~30 sec |
| CREST sampling | 30-60 min | Yes (--T option) | Generates full ensemble |
| NWChem DFT optimization | 5-15 min | Yes (separate jobs) | B3LYP/6-31G*, COSMO |
| NWChem NMR shielding | 10-30 min | Yes (separate jobs) | B3LYP/6-311+G(2d,p), COSMO |

### Total cost for 20-conformer ensemble

| Approach | Conformer Gen | DFT Opt (20×) | NMR (20×) | Total |
|----------|--------------|---------------|-----------|-------|
| **RDKit path** | ~30 sec | 100-300 min | 200-600 min | **5-15 hours** |
| **CREST path** | 30-60 min | 100-300 min | 200-600 min | **6-16 hours** |

**Parallelization strategy:**
- Conformer generation: Use all cores (RDKit numThreads=0, CREST --T)
- DFT jobs: Launch as separate Huey tasks (existing job queue)
- Each conformer DFT: NWChem can use multiple cores

**Confidence:** HIGH - Based on typical DFT calculation times and existing NWChem performance

## Integration with Existing Stack

### Changes to Existing Components

**1. Job Queue (Huey) - MINOR CHANGE**

Current: Single NWChem job per submission
New: Multiple NWChem jobs per submission (one per conformer)

```python
# Before (v1.x): Single-conformer
@huey.task()
def calculate_nmr(job_id, mol):
    geometry = optimize_geometry(mol)
    shifts = calculate_shifts(geometry)
    return shifts

# After (v2.0): Ensemble mode
@huey.task()
def calculate_nmr_ensemble(job_id, mol, method='rdkit'):
    conformers = generate_conformers(mol, method)

    # Parallel DFT optimization (existing pattern)
    opt_tasks = [optimize_geometry_task(job_id, conf_id, conf)
                 for conf_id, conf in enumerate(conformers)]

    # Filter by energy after DFT
    optimized = [t.get() for t in opt_tasks]  # Wait for all
    filtered = filter_post_dft(optimized)

    # Parallel NMR calculations
    nmr_tasks = [calculate_shifts_task(job_id, conf_id, conf)
                 for conf_id, conf in enumerate(filtered)]

    # Boltzmann averaging
    results = [t.get() for t in nmr_tasks]
    averaged = boltzmann_average_shifts(results)

    return averaged
```

**2. NWChem Integration - NO CHANGE**

Existing NWChem workflow works as-is:
- Input: XYZ coordinates
- Output: Optimized geometry + energy OR shielding tensors
- Already handles single structures
- Just call it N times for N conformers

**3. API Endpoints - MINOR ADDITION**

Add optional parameter to existing endpoint:

```python
class NMRJobRequest(BaseModel):
    smiles: str
    solvent: str = 'chcl3'
    conformers: bool = False  # NEW: Enable ensemble mode
    conformer_method: str = 'auto'  # 'auto', 'rdkit', 'crest'

@app.post("/api/jobs")
async def create_nmr_job(request: NMRJobRequest):
    if request.conformers:
        task = calculate_nmr_ensemble.schedule((job_id, mol, request.conformer_method))
    else:
        task = calculate_nmr.schedule((job_id, mol))  # Existing v1.x behavior
```

**Confidence:** HIGH - Architecture aligns with existing patterns, minimal disruption

## What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| **MMFF energies for Boltzmann weighting** | Poor DFT correlation (ρ ~ -0.1 to -0.45), can anticorrelate | DFT energies from NWChem optimization |
| **xTB energies for final weighting** | Good for filtering (ρ ~ 0.4) but not final weighting | DFT energies from NWChem |
| **CENSO workflow** | Excellent tool but requires ORCA/Turbomole | Custom filtering with NWChem |
| **OpenBabel for XYZ parsing** | Unnecessary dependency, simple format | Python stdlib (read/split) |
| **xtb-python package** | Not needed, call CREST binary | subprocess.run(['crest', ...]) |
| **scipy.stats.boltzmann** | Wrong distribution (PMF for discrete), not conformer weighting | Custom numpy implementation |
| **RDKit ETKDG for solution NMR** | Crystal structure bias from CSD | RDKit KDG (unbiased) |
| **NWChem MD for conformers** | No search algorithms, AMBER/CHARMM not suitable | CREST or RDKit |

**Confidence:** HIGH - Based on literature comparison and technical requirements

## Version Compatibility Matrix

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| **RDKit** | 2015.09.1 | >=2025.09.3 | KDG since 2015, allene/cumulene support 2025.09.1 |
| **scipy** | 1.17.0 | >=1.17.0 | Already in pyproject.toml |
| **numpy** | (any) | (via scipy) | Transitive dependency, no explicit requirement |
| **CREST** | 3.0.2 | Latest continuous | 3.0+ has integrated tblite |
| **xTB** | 6.7.1 | 6.7.1 | Required by CREST for some features |
| **NWChem** | 7.2.0 | (existing) | No changes needed |
| **Python** | 3.11 | >=3.11 | Project already requires 3.11+ |

**Python version:** Project pyproject.toml specifies `requires-python = ">=3.11"` - confirmed compatible with all components.

**Confidence:** HIGH - All versions verified from official sources

## Installation Instructions

### Development Environment

**Step 1: Verify existing dependencies**

```bash
# RDKit and scipy already in pyproject.toml
pip install -e .
```

**Step 2: Install CREST/xTB (optional)**

```bash
# Option A: Conda (recommended)
conda install -c conda-forge crest xtb

# Option B: Binary download
# Download from https://github.com/crest-lab/crest/releases
# Download from https://github.com/grimme-lab/xtb/releases
# Extract and add to PATH

# Verify installation
which crest
which xtb
crest --version
xtb --version
```

**Step 3: Set environment variable (if using OpenBLAS)**

```bash
# Avoid OpenBLAS warnings
export OPENBLAS_NUM_THREADS=1

# Add to .env or shell profile
echo "export OPENBLAS_NUM_THREADS=1" >> ~/.bashrc
```

### Production Deployment

**Dockerfile example:**

```dockerfile
FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3.11 \
    nwchem \
    conda \
    && rm -rf /var/lib/apt/lists/*

# Install CREST/xTB via conda
RUN conda install -c conda-forge crest xtb

# Install Python dependencies
COPY pyproject.toml .
RUN pip install -e .

# Set environment
ENV OPENBLAS_NUM_THREADS=1
ENV PATH="/opt/conda/bin:$PATH"
```

**Health check:**

```python
def check_conformer_capabilities():
    """Check which conformer generation methods are available."""
    import shutil
    from rdkit import Chem
    from rdkit.Chem import AllChem

    capabilities = {
        'rdkit_kdg': hasattr(AllChem, 'KDG'),
        'crest': shutil.which('crest') is not None,
        'xtb': shutil.which('xtb') is not None,
    }

    return capabilities

# Example output:
# {'rdkit_kdg': True, 'crest': True, 'xtb': True}  # Full capabilities
# {'rdkit_kdg': True, 'crest': False, 'xtb': False}  # Baseline only
```

**Confidence:** HIGH - Installation procedures verified from official documentation

## Sources

### Official Documentation
- [RDKit rdDistGeom module](https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html)
- [RDKit AllChem module](https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html)
- [CREST Documentation](https://crest-lab.github.io/crest-docs/)
- [CREST Command Line Keywords](https://crest-lab.github.io/crest-docs/page/documentation/keywords.html)
- [CREST Installation Guide](https://crest-lab.github.io/crest-docs/page/installation/install_basic.html)
- [xTB Setup and Installation](https://xtb-docs.readthedocs.io/en/latest/setup.html)

### Repositories and Releases
- [CREST GitHub Repository](https://github.com/crest-lab/crest)
- [CREST Releases](https://github.com/crest-lab/crest/releases)
- [xTB GitHub Repository](https://github.com/grimme-lab/xtb)
- [xTB Releases](https://github.com/grimme-lab/xtb/releases)
- [RDKit Releases](https://github.com/rdkit/rdkit/releases)
- [RDKit 2025.09.1 Release Discussion](https://github.com/rdkit/rdkit/discussions/8831)

### Scientific Literature
- Riniker, S. & Landrum, G. A. (2015) "Better Informed Distance Geometry: Using What We Know To Improve Conformation Generation" *J. Chem. Inf. Model.* 55(12):2562-2574 [DOI: 10.1021/acs.jcim.5b00654](https://pubs.acs.org/doi/10.1021/acs.jcim.5b00654)
- Pracht, P. et al. (2024) "CREST—A program for the exploration of low-energy molecular chemical space" *J. Chem. Phys.* 160:114110 [DOI: 10.1063/5.0197592](https://pubs.aip.org/aip/jcp/article/160/11/114110/3278084/)
- Bannwarth, C. et al. (2019) "GFN2-xTB—An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method" *J. Chem. Theory Comput.* 15:1652-1671

### Community Resources
- [RDKit Discussion #8226: Best practices - Conformer generation](https://github.com/rdkit/rdkit/discussions/8226)
- [Oxford Protein Informatics Group: Advances in Conformer Generation](https://www.blopig.com/blog/2016/06/advances-in-conformer-generation-etkdg-and-etdg/)
- [CREST File Formats Documentation](https://crest-lab.github.io/crest-docs/page/documentation/coords.html)
- [Evaluating Boltzmann Weighting Error](https://corinwagen.github.io/public/blog/20221228_boltzmann_error.html)

### Conda Packages
- [CREST on conda-forge](https://anaconda.org/conda-forge/crest)
- [xTB on conda-forge](https://anaconda.org/conda-forge/xtb)

---

*Stack research for v2.0: Conformational Sampling with Boltzmann-weighted NMR*
*Milestone: Adding ensemble NMR predictions to existing web service*
*Researched: 2026-01-26*
*Confidence: HIGH*
