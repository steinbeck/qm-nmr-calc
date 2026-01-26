"""Canonical atom ordering for consistent indexing across conformer lifecycle.

This module provides functions to establish canonical atom ordering that remains
consistent across RDKit conformer generation, NWChem optimization, and NMR parsing.

Key insight: Atom ordering is critical for correctness (Pitfall 4). When RDKit
generates conformers and NWChem optimizes them, atom indices must map back to
the original SMILES atoms for correct shift assignment.

Functions:
    canonical_atom_order: Get deterministic canonical ranking for a SMILES
    map_nwchem_to_canonical: Map 1-based NWChem indices to canonical ranks
    get_atom_count: Count atoms by type for validation
"""
from rdkit import Chem


def _smiles_to_mol_with_h(smiles: str) -> Chem.Mol:
    """Parse SMILES and add explicit hydrogens.

    Helper function to avoid duplication across module functions.

    Args:
        smiles: SMILES string representing the molecule

    Returns:
        RDKit Mol object with explicit hydrogens

    Raises:
        ValueError: If SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.AddHs(mol)


def canonical_atom_order(smiles: str) -> list[int]:
    """Get canonical atom ordering for a SMILES string.

    Creates an RDKit mol with explicit hydrogens, then returns canonical
    atom rankings using RDKit's built-in canonicalization algorithm.
    This ordering is deterministic - same SMILES always produces same ranking.

    Args:
        smiles: SMILES string representing the molecule

    Returns:
        List where position i contains the canonical rank of atom i.
        Ranks are 0-based integers from 0 to (num_atoms - 1).

    Raises:
        ValueError: If SMILES string is invalid

    Example:
        >>> order = canonical_atom_order("CCO")  # ethanol
        >>> len(order)  # 2C + 1O + 6H = 9 atoms
        9
        >>> set(order)  # All ranks from 0 to 8
        {0, 1, 2, 3, 4, 5, 6, 7, 8}
    """
    # Parse SMILES and add explicit hydrogens
    mol = _smiles_to_mol_with_h(smiles)

    # Get canonical atom ranking
    # This uses RDKit's canonical ordering algorithm, which is deterministic
    canonical_ranks = list(Chem.CanonicalRankAtoms(mol))

    return canonical_ranks


def map_nwchem_to_canonical(smiles: str) -> dict[int, int]:
    """Map 1-based NWChem atom indices to canonical ranks.

    NWChem uses 1-based indexing in its output files. This function provides
    a mapping from those 1-based indices to the canonical ranks.

    The mapping works because RDKit's atom ordering after AddHs is deterministic
    for a given SMILES, and NWChem uses the atom order from the XYZ file
    (which comes from RDKit).

    Args:
        smiles: SMILES string representing the molecule

    Returns:
        Dictionary mapping 1-based NWChem index to 0-based canonical rank.
        Keys are integers 1 through num_atoms (NWChem indices).
        Values are integers 0 through (num_atoms - 1) (canonical ranks).

    Raises:
        ValueError: If SMILES string is invalid

    Example:
        >>> mapping = map_nwchem_to_canonical("C")  # methane
        >>> len(mapping)  # 5 atoms
        5
        >>> mapping[1]  # NWChem atom 1 has canonical rank...
        0  # (example - actual rank depends on canonicalization)
    """
    # Get canonical ordering
    canonical_ranks = canonical_atom_order(smiles)

    # Create 1-based mapping for NWChem
    # NWChem uses 1-based indices, so we map from 1..N to canonical ranks
    mapping = {
        nwchem_idx: canonical_ranks[nwchem_idx - 1]
        for nwchem_idx in range(1, len(canonical_ranks) + 1)
    }

    return mapping


def get_atom_count(smiles: str) -> tuple[int, int, int]:
    """Get atom counts for a SMILES string.

    Useful for validation that atom counts match across the conformer lifecycle
    (generation -> optimization -> NMR parsing).

    Args:
        smiles: SMILES string representing the molecule

    Returns:
        Tuple of (total_atoms, num_hydrogens, num_heavy_atoms)
        - total_atoms: Total number of atoms including explicit hydrogens
        - num_hydrogens: Number of hydrogen atoms
        - num_heavy_atoms: Number of non-hydrogen atoms

    Raises:
        ValueError: If SMILES string is invalid

    Example:
        >>> get_atom_count("CCO")  # ethanol
        (9, 6, 3)  # 9 total: 2C + 1O + 6H, 6 H, 3 heavy
    """
    # Parse SMILES and add explicit hydrogens
    mol = _smiles_to_mol_with_h(smiles)

    # Count atoms
    total_atoms = mol.GetNumAtoms()

    # Count hydrogens
    num_h = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)

    # Heavy atoms (non-hydrogen)
    num_heavy = total_atoms - num_h

    return total_atoms, num_h, num_heavy
