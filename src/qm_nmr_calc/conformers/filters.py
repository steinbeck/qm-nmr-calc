"""Conformer filtering functions for deduplication and energy-based selection."""

from rdkit import Chem
from rdkit.Chem import rdMolAlign


def deduplicate_by_rmsd(mol: Chem.Mol, rmsd_threshold: float = 0.5) -> list[int]:
    """Remove duplicate conformers using symmetry-aware RMSD comparison.

    Uses a greedy algorithm: keep the first conformer, then for each subsequent
    conformer, compare it to all kept conformers. If RMSD to any kept conformer
    is below threshold, mark as duplicate and skip.

    Args:
        mol: RDKit Mol object with one or more conformers
        rmsd_threshold: RMSD cutoff in Angstroms for considering conformers identical

    Returns:
        List of conformer IDs to keep (subset of original conformer IDs)

    Example:
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> AllChem.EmbedMultipleConfs(mol, numConfs=10)
        >>> kept_ids = deduplicate_by_rmsd(mol, rmsd_threshold=0.5)
        >>> # Create new mol with only kept conformers
        >>> for conf_id in range(mol.GetNumConformers() - 1, -1, -1):
        ...     if conf_id not in kept_ids:
        ...         mol.RemoveConformer(conf_id)
    """
    # Get all conformer IDs
    conf_ids = [conf.GetId() for conf in mol.GetConformers()]

    # Single or no conformers - nothing to deduplicate
    if len(conf_ids) <= 1:
        return conf_ids

    # Greedy deduplication
    kept_ids = [conf_ids[0]]  # Always keep first

    for conf_id in conf_ids[1:]:
        is_duplicate = False

        # Compare to all kept conformers
        for kept_id in kept_ids:
            # Use GetBestRMS for symmetry-aware RMSD
            rmsd = rdMolAlign.GetBestRMS(mol, mol, prbId=conf_id, refId=kept_id)

            if rmsd < rmsd_threshold:
                is_duplicate = True
                break

        if not is_duplicate:
            kept_ids.append(conf_id)

    return kept_ids


def filter_by_energy_window(
    conf_ids: list[int], energies: list[float], window_kcal: float = 6.0
) -> tuple[list[int], list[float]]:
    """Filter conformers by energy window relative to minimum.

    Keeps only conformers within window_kcal of the lowest energy conformer.
    This is intentionally decoupled from RDKit Mol objects to make it reusable
    for both MMFF and DFT energy filtering.

    Args:
        conf_ids: List of conformer IDs
        energies: List of energies (parallel to conf_ids, same length)
        window_kcal: Energy window in kcal/mol above minimum energy

    Returns:
        Tuple of (filtered_conf_ids, filtered_energies) preserving input order

    Example:
        >>> conf_ids = [0, 1, 2, 3]
        >>> energies = [10.0, 12.0, 14.0, 20.0]
        >>> kept_ids, kept_energies = filter_by_energy_window(
        ...     conf_ids, energies, window_kcal=6.0
        ... )
        >>> kept_ids
        [0, 1, 2]
        >>> kept_energies
        [10.0, 12.0, 14.0]
    """
    # Handle empty input
    if not conf_ids or not energies:
        return ([], [])

    # Find minimum energy
    min_energy = min(energies)

    # Filter conformers within window
    filtered_ids = []
    filtered_energies = []

    for conf_id, energy in zip(conf_ids, energies):
        if (energy - min_energy) <= window_kcal:
            filtered_ids.append(conf_id)
            filtered_energies.append(energy)

    return (filtered_ids, filtered_energies)
