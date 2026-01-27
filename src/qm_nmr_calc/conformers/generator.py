"""
Conformer generation using RDKit distance geometry (KDG) and MMFF optimization.

This module provides the core conformer generation engine for the ensemble pipeline:
- KDG (distance geometry) embedding without crystal structure bias
- MMFF force field optimization with energy calculation
- Adaptive conformer count based on molecular flexibility
"""

from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdMolDescriptors


def calculate_num_conformers(smiles: str, max_conformers: Optional[int] = None) -> int:
    """
    Calculate adaptive number of conformers based on molecular flexibility.

    Strategy:
    - Rigid molecules (≤8 rotatable bonds): 50 conformers
    - Flexible molecules (>8 rotatable bonds): 200 conformers
    - Override with max_conformers if provided

    Args:
        smiles: SMILES string of molecule
        max_conformers: Optional override for conformer count

    Returns:
        Number of conformers to generate

    Raises:
        ValueError: If SMILES is invalid
    """
    # If override provided, use it directly
    if max_conformers is not None:
        return max_conformers

    # Parse SMILES without hydrogens to count rotatable bonds
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Count rotatable bonds
    num_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Return adaptive count
    return 200 if num_rotatable > 8 else 50


def generate_conformers_kdg(
    smiles: str, num_confs: int, random_seed: int = 0xF00D
) -> Chem.Mol:
    """
    Generate conformers using distance geometry with KDG parameters.

    Uses RDKit's KDG (Knowledge-based Distance Geometry) without ETKDG's
    crystal structure bias, appropriate for solution-phase NMR predictions.

    Configuration:
    - No pre-optimization pruning (pruneRmsThresh=-1.0)
    - All available threads (numThreads=0)
    - Random coordinate fallback for difficult embeddings

    Args:
        smiles: SMILES string of molecule
        num_confs: Number of conformers to generate
        random_seed: Random seed for reproducibility (default: 0xF00D)

    Returns:
        RDKit Mol object with embedded 3D conformers (including hydrogens)

    Raises:
        ValueError: If SMILES is invalid
        RuntimeError: If conformer generation fails even with random coords
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Add hydrogens before embedding
    mol = Chem.AddHs(mol)

    # Configure KDG parameters
    params = rdDistGeom.KDG()
    params.randomSeed = random_seed
    params.numThreads = 0  # Use all threads
    params.pruneRmsThresh = -1.0  # No pre-optimization pruning

    # Attempt embedding
    conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    # If no conformers generated, retry with random coordinates
    if len(conf_ids) == 0:
        params.useRandomCoords = True
        conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

        # If still no conformers, raise error
        if len(conf_ids) == 0:
            raise RuntimeError(
                f"Failed to generate conformers for {smiles} even with random coordinates"
            )

    return mol


def optimize_conformers_mmff(
    mol: Chem.Mol, max_iters: int = 200
) -> list[tuple[int, float]]:
    """
    Optimize conformers using MMFF force field and return energies.

    MMFF (Merck Molecular Force Field) provides energies in kcal/mol.
    Optimization uses all available threads for parallel processing.

    Args:
        mol: RDKit Mol object with conformers
        max_iters: Maximum optimization iterations per conformer

    Returns:
        List of (not_converged, energy) tuples, one per conformer.
        - not_converged: 0 if converged, 1 if not converged
        - energy: MMFF energy in kcal/mol

    Raises:
        ValueError: If molecule lacks MMFF parameters
    """
    # Check if molecule has MMFF parameters
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError(
            "Molecule does not have all required MMFF parameters. "
            "This typically occurs with unusual atom types or valences."
        )

    # Optimize all conformers
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=max_iters)

    return results
