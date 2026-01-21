"""Geometry handling for NWChem calculations.

Provides SMILES-to-3D conversion using RDKit's ETKDGv3 method
and loading of pre-optimized geometries from XYZ/SDF files.
"""

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds


def smiles_to_xyz(smiles: str, random_seed: int = 0xF00D) -> tuple[Chem.Mol, str]:
    """Convert SMILES to 3D coordinates using ETKDGv3 and UFF optimization.

    Args:
        smiles: SMILES string representing the molecule
        random_seed: Random seed for reproducible conformer generation

    Returns:
        Tuple of (RDKit Mol object with 3D coordinates, XYZ coordinate block string)

    Raises:
        ValueError: If SMILES string is invalid
        RuntimeError: If 3D embedding fails
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Add hydrogens (required for realistic 3D geometry)
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates with ETKDGv3 (best method per RESEARCH.md)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    result = AllChem.EmbedMolecule(mol, params)

    if result != 0:
        raise RuntimeError(
            f"Failed to embed molecule from SMILES '{smiles}'. "
            "The molecule may be too constrained or have unusual valences."
        )

    # Run UFF force field optimization
    AllChem.UFFOptimizeMolecule(mol)

    # Convert to XYZ block for NWChem input
    xyz_block = mol_to_xyz_block(mol)

    return mol, xyz_block


def load_geometry_file(
    filepath: Path | str, charge: int = 0
) -> tuple[Chem.Mol, str]:
    """Load molecular geometry from XYZ or SDF file.

    Args:
        filepath: Path to .xyz or .sdf file
        charge: Molecular formal charge (required for XYZ bond determination)

    Returns:
        Tuple of (RDKit Mol object with 3D coordinates, XYZ coordinate block string)

    Raises:
        ValueError: If file doesn't parse or has unsupported extension
        RuntimeError: If bond determination fails for XYZ files
    """
    filepath = Path(filepath)
    ext = filepath.suffix.lower()

    if ext == ".sdf":
        # SDF contains bond info, straightforward
        mol = Chem.MolFromMolFile(str(filepath), removeHs=False)
        if mol is None:
            raise ValueError(f"Failed to parse SDF file: {filepath}")

    elif ext == ".xyz":
        # XYZ only has coordinates, must determine bonds
        mol = Chem.MolFromXYZFile(str(filepath))
        if mol is None:
            raise ValueError(f"Failed to parse XYZ file: {filepath}")

        # Determine bonds from coordinates (requires charge!)
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=charge)
        except Exception as e:
            raise RuntimeError(
                f"Bond determination failed for {filepath}: {e}. "
                "Check that the charge parameter is correct."
            ) from e

    else:
        raise ValueError(
            f"Unsupported file format: {ext}. Supported formats: .xyz, .sdf"
        )

    # Validate geometry
    validate_geometry(mol)

    # Convert to XYZ block for NWChem input
    xyz_block = mol_to_xyz_block(mol)

    return mol, xyz_block


def mol_to_xyz_block(mol: Chem.Mol) -> str:
    """Convert RDKit mol to XYZ coordinate block for NWChem geometry directive.

    The output contains only the "Element X Y Z" lines, without the atom count
    and title lines that are in standard XYZ format.

    Args:
        mol: RDKit Mol object with 3D coordinates

    Returns:
        XYZ coordinate block string (element + 3 coordinates per line)

    Raises:
        ValueError: If molecule has no conformer or coordinates
    """
    # Validate before conversion
    validate_geometry(mol)

    # Get full XYZ block from RDKit
    full_xyz = Chem.MolToXYZBlock(mol)

    # Skip first two lines (atom count and title/comment)
    lines = full_xyz.strip().split("\n")
    if len(lines) < 3:
        raise ValueError("XYZ block too short - no atom coordinates found")

    # Return just the coordinate lines
    return "\n".join(lines[2:])


def validate_geometry(mol: Chem.Mol) -> bool:
    """Validate that a molecule has valid 3D coordinates.

    Args:
        mol: RDKit Mol object to validate

    Returns:
        True if valid

    Raises:
        ValueError: If molecule has no conformer or invalid coordinates
    """
    if mol is None:
        raise ValueError("Molecule is None")

    # Check mol has at least one conformer
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no 3D conformer")

    # Get the conformer
    conf = mol.GetConformer(0)

    # Check that coordinates are not all zeros (sign of failed embedding)
    all_zeros = True
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        if pos.x != 0.0 or pos.y != 0.0 or pos.z != 0.0:
            all_zeros = False
            break

    if all_zeros:
        raise ValueError("Molecule has invalid coordinates (all zeros)")

    return True
