# Conformational Sampling and Boltzmann-Weighted NMR Calculations for qm-nmr-calc

## Executive Summary

This document provides a comprehensive analysis of approaches for implementing conformational sampling with Boltzmann-weighted averaging for NMR chemical shift calculations in your qm-nmr-calc project. Based on extensive research of current methodologies, I recommend a **tiered approach** using CREST/xTB for conformer generation followed by DFT refinement in NWChem, with RDKit as a faster alternative for simpler molecules.

---

## 1. Why Conformational Sampling Matters

For flexible molecules, NMR chemical shifts represent a population-weighted average over all accessible conformations at the measurement temperature. The ISiCLE paper demonstrated that single-conformer calculations can produce outliers exceeding 15 ppm from experimental values, whereas Boltzmann-weighted ensembles significantly reduce prediction errors.

### Key Evidence

The ISiCLE authors tested methylcyclohexane with 80 conformers (40 axial, 40 equatorial), obtaining a relative free energy difference of 1.99 kcal/mol between conformer classes using NWChem—consistent with high-level Gaussian calculations and experimental data.

**Critical insight**: At 298 K, conformers within ~2.5 kcal/mol of the global minimum have significant Boltzmann populations. A conformer at 3 kcal/mol above the minimum contributes only ~0.6% to the weighted average, but one at 1 kcal/mol contributes ~18%.

---

## 2. Approach Comparison

### 2.1 Option A: RDKit Conformer Generation (Recommended for Simpler Molecules)

**Workflow**: RDKit ETKDG → MMFF optimization → DFT optimization in NWChem → NMR calculation

**Advantages**:
- Already integrated in your project (RDKit dependency exists)
- Fast conformer enumeration using ETKDGv3
- Good for drug-like molecules with moderate flexibility
- No additional dependencies required

**Limitations**:
- ETKDGv3 uses CSD-derived torsion preferences (biased toward crystal structures)
- For solution-phase NMR, pure distance geometry (KDG) may sample better
- May miss important conformers for highly flexible molecules
- MMFF energy ranking has poor correlation with DFT (Spearman ρ ~ -0.1 to -0.45)

**Implementation**:
```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import numpy as np

def generate_conformers_rdkit(mol, num_confs=50, energy_window=5.0, rms_thresh=0.5):
    """Generate conformer ensemble using RDKit."""
    mol = Chem.AddHs(mol)
    
    # Use KDG for solution sampling (no crystal bias)
    # or ETKDGv3() for molecules where crystal-like conformers are relevant
    params = AllChem.KDGParams()  # Pure distance geometry
    params.pruneRmsThresh = rms_thresh
    params.numThreads = 0  # Use all cores
    
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    
    # MMFF optimization and energy calculation
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant='MMFF94s')
    
    # Filter by energy window
    energies = [r[1] for r in results if r[0] == 0]  # Only converged
    min_energy = min(energies)
    
    filtered_confs = []
    for i, (converged, energy) in enumerate(results):
        if converged == 0 and (energy - min_energy) <= energy_window:
            filtered_confs.append((cids[i], energy))
    
    return mol, filtered_confs
```

### 2.2 Option B: CREST/xTB (Recommended for Accurate Results)

**Workflow**: CREST metadynamics with GFN2-xTB → DFT optimization in NWChem → NMR calculation

**Advantages**:
- State-of-the-art conformer sampling using metadynamics
- GFN2-xTB provides much better conformer ranking than force fields (Spearman ρ ~ 0.39-0.47 vs DFT)
- Comprehensive sampling even for macrocycles and complex systems
- Well-validated for NMR applications (Grimme lab's CENSO workflow)
- Solvation models available (ALPB, GBSA)

**Limitations**:
- Additional dependency (xtb, crest binaries)
- Slower than RDKit (~hours for drug-like molecules vs seconds)
- Overkill for rigid molecules

**Implementation**:
```python
import subprocess
import os
from pathlib import Path

def generate_conformers_crest(xyz_file, charge=0, solvent=None, ewin=6.0, temp=298.15):
    """Generate conformer ensemble using CREST/xTB."""
    
    cmd = [
        'crest', xyz_file,
        '--gfn2',
        '--chrg', str(charge),
        '--ewin', str(ewin),
        '--T', str(os.cpu_count() or 4),
        '-temp', str(temp)
    ]
    
    if solvent:
        cmd.extend(['--alpb', solvent])
    
    subprocess.run(cmd, check=True)
    
    # Output: crest_conformers.xyz
    return parse_crest_ensemble('crest_conformers.xyz')

def parse_crest_ensemble(ensemble_file):
    """Parse CREST multi-xyz output."""
    conformers = []
    with open(ensemble_file) as f:
        content = f.read()
    
    # Split by number of atoms line (CREST format)
    blocks = content.strip().split('\n\n')
    for i, block in enumerate(blocks):
        lines = block.strip().split('\n')
        n_atoms = int(lines[0])
        # Energy is in the comment line
        energy = float(lines[1].split()[0])  # Hartree
        coords = '\n'.join(lines[2:n_atoms+2])
        conformers.append({
            'id': i,
            'energy_hartree': energy,
            'xyz': coords
        })
    
    return conformers
```

### 2.3 Option C: NWChem Native MD (Not Recommended)

NWChem does have molecular dynamics capabilities, but:
- Classical MD module uses AMBER/CHARMM force fields (not QM)
- Car-Parrinello MD is for periodic systems with plane-wave basis
- No built-in conformer search algorithms
- Would require extensive custom scripting

**Verdict**: Use NWChem only for the final DFT optimization and NMR calculation steps.

---

## 3. Boltzmann Weighting Implementation

### 3.1 The Boltzmann Formula

For conformers at thermal equilibrium:

```
p_i = exp(-E_i / RT) / Σ exp(-E_j / RT)
```

Where:
- p_i = population of conformer i
- E_i = relative energy (typically Gibbs free energy)
- R = 8.314 J/(mol·K) = 0.001987 kcal/(mol·K)
- T = temperature (typically 298.15 K)

The Boltzmann-averaged chemical shift:

```
δ_avg = Σ p_i × δ_i
```

### 3.2 Energy Source for Weighting

**Critical decision**: Which energies to use for Boltzmann weighting?

| Level | Source | Accuracy | Notes |
|-------|--------|----------|-------|
| Low | MMFF (RDKit) | Poor | May anticorrelate with DFT |
| Medium | GFN2-xTB | Good (~2 kcal/mol MAE) | Fast, good for filtering |
| High | DFT single-point | Very good | Use same level as NMR calc |
| Highest | DFT + thermal (RRHO) | Excellent | Includes entropy corrections |

**Recommendation**: Use DFT single-point energies from the geometry optimization step for weighting. If accuracy is critical, add GFN2-xTB quasi-harmonic RRHO corrections for entropy.

### 3.3 Python Implementation

```python
import numpy as np

def boltzmann_weights(energies_kcal, temperature=298.15):
    """
    Calculate Boltzmann weights from relative energies.
    
    Args:
        energies_kcal: Array of relative energies in kcal/mol (min should be 0)
        temperature: Temperature in Kelvin
    
    Returns:
        Array of normalized Boltzmann weights (sum = 1)
    """
    R = 0.001987204  # kcal/(mol·K)
    RT = R * temperature
    
    # Ensure minimum energy is 0
    rel_energies = energies_kcal - np.min(energies_kcal)
    
    # Calculate Boltzmann factors
    boltz_factors = np.exp(-rel_energies / RT)
    
    # Normalize
    weights = boltz_factors / np.sum(boltz_factors)
    
    return weights

def boltzmann_average_shifts(shifts_per_conformer, energies_kcal, temperature=298.15):
    """
    Calculate Boltzmann-weighted average chemical shifts.
    
    Args:
        shifts_per_conformer: List of dicts, each mapping atom_idx -> shift
        energies_kcal: Array of conformer energies
        temperature: Temperature in Kelvin
    
    Returns:
        Dict mapping atom_idx -> averaged shift
    """
    weights = boltzmann_weights(energies_kcal, temperature)
    
    # Get all atom indices
    all_atoms = set()
    for shifts in shifts_per_conformer:
        all_atoms.update(shifts.keys())
    
    # Calculate weighted average for each atom
    avg_shifts = {}
    for atom_idx in all_atoms:
        weighted_sum = 0.0
        weight_sum = 0.0
        
        for i, shifts in enumerate(shifts_per_conformer):
            if atom_idx in shifts:
                weighted_sum += weights[i] * shifts[atom_idx]
                weight_sum += weights[i]
        
        if weight_sum > 0:
            avg_shifts[atom_idx] = weighted_sum / weight_sum
    
    return avg_shifts, weights
```

---

## 4. Recommended Workflow for qm-nmr-calc

### 4.1 Complete Pipeline

```
Input SMILES/MOL
       ↓
[1] 3D Embedding (RDKit)
       ↓
[2] Conformer Generation
    ├─ Simple molecules (≤3 rotatable bonds): RDKit ETKDG
    └─ Complex molecules: CREST/xTB
       ↓
[3] Energy Filtering (keep ΔE < 5-6 kcal/mol)
       ↓
[4] DFT Geometry Optimization (NWChem)
    - B3LYP/6-31G* with COSMO
    - Rerank by DFT energies
       ↓
[5] Energy Filtering (keep ΔE < 3 kcal/mol or top 90% Boltzmann)
       ↓
[6] NMR Shielding Calculation (NWChem)
    - B3LYP/6-311+G(2d,p) GIAO with COSMO
       ↓
[7] Boltzmann Averaging
    - Weight by DFT optimization energies
       ↓
[8] Convert to Chemical Shifts
    - δ = σ_TMS - σ_sample (or linear scaling)
```

### 4.2 Integration with Existing Code

Your current qm-nmr-calc already handles:
- SMILES → 3D geometry (RDKit)
- NWChem input generation
- NMR shielding calculation
- Shift conversion with scaling

**New components needed**:

1. **Conformer generation module** (`src/qm_nmr_calc/conformers.py`)
2. **Ensemble job management** (modify `tasks.py`)
3. **Boltzmann averaging** (`src/qm_nmr_calc/averaging.py`)

### 4.3 NWChem Input Template for Ensemble

```python
def generate_nwchem_conformer_opt(xyz_coords, conf_id, solvent='chcl3'):
    """Generate NWChem input for conformer optimization."""
    
    SOLVENT_EPSILON = {'chcl3': 4.81, 'dmso': 46.7}
    eps = SOLVENT_EPSILON.get(solvent, 4.81)
    
    return f"""
title "Conformer {conf_id} optimization"
start conf_{conf_id}

geometry units angstroms
{xyz_coords}
end

basis
  * library 6-31G*
end

dft
  xc b3lyp
  mult 1
end

cosmo
  dielec {eps}
end

driver
  maxiter 100
  tight
end

task dft optimize
"""

def generate_nwchem_nmr_ensemble(optimized_xyz, conf_id, solvent='chcl3'):
    """Generate NWChem input for NMR calculation."""
    
    SOLVENT_EPSILON = {'chcl3': 4.81, 'dmso': 46.7}
    eps = SOLVENT_EPSILON.get(solvent, 4.81)
    
    return f"""
title "Conformer {conf_id} NMR"
start nmr_{conf_id}

geometry units angstroms
{optimized_xyz}
end

basis
  * library 6-311+G(2d,p)
end

dft
  xc b3lyp
  mult 1
end

cosmo
  dielec {eps}
end

property
  shielding
end

task dft property
"""
```

---

## 5. Practical Considerations

### 5.1 Conformer Count Guidelines

| Molecule Type | Rotatable Bonds | Recommended Conformers | Method |
|---------------|-----------------|----------------------|--------|
| Rigid | 0-1 | 1-5 | RDKit quick |
| Semi-rigid | 2-4 | 10-30 | RDKit ETKDG |
| Flexible | 5-8 | 30-100 | CREST or extensive RDKit |
| Very flexible | >8 | 100-500 | CREST mandatory |

### 5.2 Energy Window Selection

- **Initial filtering**: 6 kcal/mol (captures 99%+ of population at 298K)
- **Post-DFT filtering**: 3 kcal/mol or keep conformers representing 95-99% cumulative Boltzmann population
- **For NMR**: Keep conformers with population > 1%

### 5.3 Computational Cost Estimation

For a typical drug-like molecule (~30 atoms):

| Step | Time (est.) | Parallelizable |
|------|-------------|----------------|
| CREST sampling | 30-60 min | Yes (xTB threads) |
| DFT opt per conformer | 5-15 min | Yes (separate jobs) |
| NMR per conformer | 10-30 min | Yes (separate jobs) |

For 20 conformers: ~6-15 hours total (with parallelization)

### 5.4 Error Analysis

Research shows that Boltzmann weighting is remarkably robust to energy errors:
- With 0.28 kcal/mol random energy errors, the weighted NMR shift was closer to the true value than the single-conformer shift 99% of the time
- Even with substantial fluctuations in major conformer weights, the averaged property remains stable

---

## 6. References and Resources

### Key Papers

1. **ISiCLE/NWChem workflow**: Yesiltepe et al. (2018) "An automated framework for NMR chemical shift calculations of small organic molecules" *J. Cheminformatics* 10:52

2. **CREST methodology**: Pracht et al. (2020) "Automated exploration of the low-energy chemical space with fast quantum chemical methods" *Phys. Chem. Chem. Phys.* 22:7169-7192

3. **GFN2-xTB**: Bannwarth et al. (2019) "GFN2-xTB—An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method" *J. Chem. Theory Comput.* 15:1652-1671

4. **DP4+ NMR validation**: Grimblat et al. (2017) "Beyond DP4: an Improved Probability for Stereochemical Assignment" *J. Org. Chem.* 80:12526-12534

5. **Force field comparison for NMR**: Cuadrado et al. (2022) "May the Force (Field) Be with You: On the Importance of Conformational Searches" *Mar. Drugs* 20:699

6. **CENSO workflow**: Grimme et al. (2021) "Efficient quantum chemical calculation of structure ensembles and free energies for non-rigid molecules" *J. Phys. Chem. A* 125:4039-4054

### Software Links

- CREST: https://github.com/crest-lab/crest
- xTB: https://github.com/grimme-lab/xtb
- CENSO: https://github.com/grimme-lab/CENSO
- RDKit: https://www.rdkit.org

---

## 7. Implementation Roadmap

### Phase 1: Minimal Viable Feature
1. Add RDKit conformer generation
2. Implement sequential NWChem calculations for each conformer
3. Add Boltzmann averaging module
4. Update API to support conformer sampling option

### Phase 2: Improved Accuracy
1. Integrate CREST/xTB as optional backend
2. Add conformer filtering by Boltzmann population
3. Implement parallel conformer processing
4. Add thermal corrections (RRHO)

### Phase 3: Production Quality
1. Add adaptive conformer count based on molecule flexibility
2. Implement caching of conformer ensembles
3. Add DP4+ probability calculation for structure validation
4. Performance optimization and benchmarking

---

## 8. Conclusion

For qm-nmr-calc, I recommend:

1. **Start with RDKit** for conformer generation (already a dependency)
2. **Use KDG (not ETKDG)** for solution-phase NMR to avoid crystal bias
3. **Add CREST support** as an optional high-accuracy mode
4. **Weight by DFT energies** from the optimization step
5. **Use 3-5 kcal/mol energy window** for initial filtering

This approach balances accuracy with computational cost and leverages your existing NWChem infrastructure.
