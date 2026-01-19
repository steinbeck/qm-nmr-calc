"""Molecule validation using RDKit."""

import io
import sys
from typing import Optional

from rdkit import Chem, RDLogger


def validate_smiles(smiles: str) -> tuple[Optional[Chem.Mol], Optional[str]]:
    """Validate SMILES and return molecule or error message.

    Args:
        smiles: SMILES string to validate

    Returns:
        (mol, None) on success
        (None, error_message) on failure
    """
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


def validate_mol_file(
    content: bytes, filename: str
) -> tuple[Optional[Chem.Mol], Optional[str]]:
    """Validate MOL/SDF file content.

    Args:
        content: Raw bytes of file content
        filename: Original filename (used for error messages)

    Returns:
        (mol, None) on success with single molecule
        (None, error_message) on failure or multiple molecules
    """
    content_str = content.decode("utf-8", errors="replace")

    # Check if SDF (contains $$$$ delimiter)
    if "$$$$" in content_str:
        # SDF file - check for multiple molecules
        suppl = Chem.SDMolSupplier()
        suppl.SetData(content_str)

        mols = [m for m in suppl if m is not None]

        if len(mols) == 0:
            return None, "No valid molecules found in SDF file"
        if len(mols) > 1:
            return (
                None,
                f"SDF file contains {len(mols)} molecules. Only single-molecule files are accepted.",
            )

        return mols[0], None
    else:
        # MOL file
        mol = Chem.MolFromMolBlock(content_str)
        if mol is None:
            return None, "Invalid MOL file format"
        return mol, None
