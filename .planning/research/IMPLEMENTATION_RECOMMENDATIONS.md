# Implementation Recommendations: Conformer Pre-Selection

**Project:** qm-nmr-calc conformer pipeline enhancement
**Date:** 2026-01-30

---

## Recommended Approach

Based on research findings, implement a **two-stage pre-selection**:

1. **RMSD Clustering** - Reduce redundancy (RDKit-only, immediate)
2. **xTB Energy Ranking** - Better energy ordering (optional, recommended)

This should reduce ~40 conformers to ~8 while maintaining coverage of conformational space.

---

## Option A: RDKit-Only (Quick Implementation)

### Changes to `conformers/filters.py`

Add Butina clustering function:

```python
from rdkit.Chem import rdMolAlign
from rdkit.ML.Cluster import Butina

def cluster_conformers_by_rmsd(
    mol: Chem.Mol,
    rmsd_threshold: float = 1.5,
) -> list[list[int]]:
    """Cluster conformers using Butina algorithm with symmetry-aware RMSD.

    Args:
        mol: RDKit Mol with multiple conformers
        rmsd_threshold: RMSD threshold in Angstroms for clustering

    Returns:
        List of clusters, each cluster is a list of conformer IDs
        Clusters are sorted by size (largest first)
    """
    conf_ids = [conf.GetId() for conf in mol.GetConformers()]
    n_confs = len(conf_ids)

    if n_confs <= 1:
        return [[conf_ids[0]]] if conf_ids else []

    # Build distance matrix (lower triangle, required by Butina)
    dists = []
    for i in range(1, n_confs):
        for j in range(i):
            rmsd = rdMolAlign.GetBestRMS(
                mol, mol,
                prbId=conf_ids[i],
                refId=conf_ids[j]
            )
            dists.append(rmsd)

    # Butina clustering
    clusters = Butina.ClusterData(
        dists,
        n_confs,
        distThresh=rmsd_threshold,
        isDistData=True
    )

    # Convert indices back to conformer IDs
    return [[conf_ids[idx] for idx in cluster] for cluster in clusters]


def select_cluster_representatives(
    clusters: list[list[int]],
    energies: dict[int, float],
    max_representatives: int = 8,
) -> list[int]:
    """Select lowest-energy conformer from each cluster.

    Args:
        clusters: List of clusters from cluster_conformers_by_rmsd
        energies: Dict mapping conformer ID to energy
        max_representatives: Maximum number of representatives to return

    Returns:
        List of conformer IDs (one per cluster, lowest energy)
    """
    representatives = []

    for cluster in clusters:
        if len(representatives) >= max_representatives:
            break

        # Find lowest energy conformer in this cluster
        cluster_with_energy = [(cid, energies.get(cid, float('inf'))) for cid in cluster]
        cluster_with_energy.sort(key=lambda x: x[1])
        representatives.append(cluster_with_energy[0][0])

    return representatives
```

### Changes to `conformers/pipeline.py`

After MMFF optimization and before DFT:

```python
from .filters import cluster_conformers_by_rmsd, select_cluster_representatives

# After MMFF optimization, before creating ConformerEnsemble:

# Step 5.1: Cluster by RMSD
clusters = cluster_conformers_by_rmsd(mol, rmsd_threshold=1.5)

# Step 5.2: Build energy lookup
energy_lookup = {conf_id: energy for conf_id, energy in zip(conf_ids, energies)}

# Step 5.3: Select representatives (one per cluster, lowest MMFF energy)
target_conformers = min(8, len(clusters))
representative_ids = select_cluster_representatives(
    clusters,
    energy_lookup,
    max_representatives=target_conformers
)

# Step 5.4: Filter to only representatives
final_conf_ids = representative_ids
final_energies = [energy_lookup[cid] for cid in final_conf_ids]
```

### Expected Results

| Metric | Before | After |
|--------|--------|-------|
| Conformers to DFT | ~40 | ~8 |
| DFT time | ~10-20 hours | ~2-4 hours |
| Coverage quality | Full but redundant | Diverse representatives |

---

## Option B: xTB Integration (Better Accuracy)

### Prerequisites

Install xTB binary:
```bash
# Via conda (recommended)
conda install -c conda-forge xtb

# Or download from GitHub releases
# https://github.com/grimme-lab/xtb/releases
```

### New module: `conformers/xtb_ranking.py`

```python
"""xTB-based conformer energy ranking for pre-DFT selection."""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import Optional
import shutil

def detect_xtb_available() -> bool:
    """Check if xTB binary is available in PATH."""
    return shutil.which("xtb") is not None


def calculate_xtb_energy(
    xyz_content: str,
    charge: int = 0,
    solvent: Optional[str] = None,
) -> float:
    """Calculate GFN2-xTB energy for a single conformer.

    Args:
        xyz_content: XYZ file content as string
        charge: Molecular charge
        solvent: ALPB solvent name (e.g., "chcl3", "dmso") or None

    Returns:
        Energy in Hartree

    Raises:
        RuntimeError: If xTB calculation fails
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        xyz_path = Path(tmpdir) / "input.xyz"
        xyz_path.write_text(xyz_content)

        cmd = [
            "xtb",
            str(xyz_path),
            "--gfn2",
            "--chrg", str(charge),
            "--sp",  # Single point only (fast)
        ]

        if solvent:
            # Map common solvent names to xTB ALPB
            solvent_map = {
                "chcl3": "chcl3",
                "chloroform": "chcl3",
                "dmso": "dmso",
            }
            alpb_solvent = solvent_map.get(solvent.lower(), solvent)
            cmd.extend(["--alpb", alpb_solvent])

        result = subprocess.run(
            cmd,
            cwd=tmpdir,
            capture_output=True,
            text=True,
            timeout=60,  # 1 minute max per conformer
        )

        if result.returncode != 0:
            raise RuntimeError(f"xTB failed: {result.stderr}")

        # Parse energy from output
        # Look for "TOTAL ENERGY" line
        for line in result.stdout.split("\n"):
            if "TOTAL ENERGY" in line:
                # Format: "          | TOTAL ENERGY              -XX.XXXXXX Eh   |"
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "TOTAL" and i + 2 < len(parts):
                        energy_str = parts[i + 2]
                        return float(energy_str)

        raise RuntimeError("Could not parse xTB energy from output")


def rank_conformers_by_xtb(
    mol,  # RDKit Mol with conformers
    conf_ids: list[int],
    charge: int = 0,
    solvent: Optional[str] = None,
) -> dict[int, float]:
    """Calculate xTB energies for multiple conformers.

    Args:
        mol: RDKit Mol with conformers
        conf_ids: List of conformer IDs to rank
        charge: Molecular charge
        solvent: Solvent name or None

    Returns:
        Dict mapping conformer ID to xTB energy in kcal/mol (relative)
    """
    from rdkit.Chem import rdmolfiles

    energies_hartree = {}

    for conf_id in conf_ids:
        # Generate XYZ content
        xyz_content = rdmolfiles.MolToXYZBlock(mol, confId=conf_id)

        try:
            energy = calculate_xtb_energy(xyz_content, charge, solvent)
            energies_hartree[conf_id] = energy
        except Exception as e:
            # Skip failed conformers, log warning
            print(f"Warning: xTB failed for conformer {conf_id}: {e}")
            continue

    if not energies_hartree:
        raise RuntimeError("All xTB calculations failed")

    # Convert to kcal/mol relative to minimum
    hartree_to_kcal = 627.509
    min_energy = min(energies_hartree.values())

    return {
        conf_id: (energy - min_energy) * hartree_to_kcal
        for conf_id, energy in energies_hartree.items()
    }
```

### Modified pipeline with xTB

```python
from .xtb_ranking import detect_xtb_available, rank_conformers_by_xtb

# In generate_conformer_ensemble():

# After MMFF opt and energy window filter:

if detect_xtb_available():
    # Use xTB for better ranking
    xtb_energies = rank_conformers_by_xtb(
        mol,
        deduped_conf_ids,
        charge=charge,
        solvent=solvent,
    )
    ranking_energies = xtb_energies
    ranking_method = "xtb_gfn2"
else:
    # Fall back to MMFF
    ranking_energies = {cid: e for cid, e in zip(deduped_conf_ids, deduped_energies)}
    ranking_method = "mmff94s"

# Cluster by RMSD
clusters = cluster_conformers_by_rmsd(mol, rmsd_threshold=1.5)

# Select representatives using the best available ranking
representatives = select_cluster_representatives(
    clusters,
    ranking_energies,
    max_representatives=8,
)
```

---

## Testing Strategy

### Unit Tests

```python
def test_rmsd_clustering_reduces_conformers():
    """Verify RMSD clustering produces fewer conformers."""
    mol = Chem.MolFromSmiles("CCCCCCC")  # Heptane
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, numConfs=50)
    AllChem.MMFFOptimizeMoleculeConfs(mol)

    clusters = cluster_conformers_by_rmsd(mol, rmsd_threshold=1.5)

    assert len(clusters) < 50  # Should reduce count
    assert len(clusters) >= 3  # But not too aggressively


def test_cluster_representatives_selects_lowest_energy():
    """Verify lowest energy conformer selected from each cluster."""
    clusters = [[0, 1, 2], [3, 4], [5]]
    energies = {0: 10.0, 1: 8.0, 2: 12.0, 3: 5.0, 4: 6.0, 5: 15.0}

    reps = select_cluster_representatives(clusters, energies, max_representatives=3)

    assert reps == [1, 3, 5]  # Lowest from each cluster
```

### Integration Test

```python
def test_pipeline_reduces_to_target_conformers():
    """Verify full pipeline reduces to ~8 conformers."""
    ensemble = generate_conformer_ensemble(
        smiles="CCCCCCCCCC",  # Decane - very flexible
        job_id="test_clustering",
        target_conformers=8,
    )

    assert len(ensemble.conformers) <= 10
    assert len(ensemble.conformers) >= 4  # Should have diversity
```

---

## Configuration Parameters

Add to calculation options:

```python
class ConformerOptions(BaseModel):
    # Existing
    max_conformers: int = 50
    energy_window_kcal: float = 6.0
    rmsd_threshold: float = 0.5  # For deduplication

    # New
    clustering_rmsd_threshold: float = 1.5  # For diversity clustering
    target_conformers_for_dft: int = 8
    use_xtb_ranking: bool = True  # Auto-detect xTB
```

---

## Rollout Plan

### Phase 1: RDKit Clustering (1-2 days)
- Add clustering functions to `filters.py`
- Update pipeline to use clustering
- Add tests
- Validate on test molecules

### Phase 2: xTB Integration (2-3 days)
- Add `xtb_ranking.py` module
- Add xTB detection and fallback
- Update pipeline to prefer xTB
- Add xTB-specific tests

### Phase 3: Tuning (1 week)
- Run on representative molecules
- Tune RMSD threshold
- Tune target conformer count
- Document performance vs accuracy tradeoffs

---

## Expected Outcomes

| Molecule Type | Current DFT Time | New DFT Time | Accuracy Impact |
|--------------|------------------|--------------|-----------------|
| Rigid | ~30 min | ~15 min | None |
| Semi-rigid (5 confs) | ~1.5 hours | ~1.5 hours | None |
| Flexible (40 confs) | ~10 hours | ~2 hours | Minimal (<0.1 ppm) |
| Very flexible (100+ confs) | ~25+ hours | ~3-4 hours | Small (<0.2 ppm) |

The accuracy impact is expected to be minimal because:
1. RMSD clustering preserves conformational diversity
2. xTB ranking correlates reasonably with DFT
3. Boltzmann averaging is robust to conformer selection
