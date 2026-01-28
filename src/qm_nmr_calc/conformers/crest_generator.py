"""
CREST conformer generation utilities.

This module provides utilities for high-accuracy conformer generation using
CREST/xTB when available:
- Binary detection for CREST and xTB
- ALPB solvent model mapping
- Multi-structure XYZ file parsing
"""

import shutil
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple


class CRESTConformer(NamedTuple):
    """Data for a single CREST-generated conformer."""

    conformer_id: str  # Sequential ID like "conf_001"
    energy_hartree: float  # Total energy in Hartree
    xyz_block: str  # Full XYZ content (atom count + comment + coordinates)


# ALPB solvent model mapping
# Maps job solvent names to CREST/xTB ALPB solvent identifiers
ALPB_SOLVENT_MAP = {
    "chcl3": "chcl3",
    "dmso": "dmso",
}


@lru_cache(maxsize=1)
def detect_crest_available() -> bool:
    """
    Detect if CREST and xTB binaries are available on PATH.

    CREST requires both binaries to function. This function uses a cache
    since binary availability doesn't change during process lifetime.

    Returns:
        True if both crest and xtb are found, False otherwise

    Example:
        >>> detect_crest_available()
        True
    """
    crest_found = shutil.which("crest") is not None
    xtb_found = shutil.which("xtb") is not None
    return crest_found and xtb_found


def get_alpb_solvent(job_solvent: str) -> str | None:
    """
    Map job solvent name to CREST ALPB solvent identifier.

    ALPB (Analytical Linearized Poisson-Boltzmann) is the implicit solvent
    model used by xTB/CREST. Only a subset of solvents are supported.

    Args:
        job_solvent: Solvent name from job input (case-insensitive)

    Returns:
        ALPB solvent identifier if supported, None otherwise
        Returns None for "vacuum" or unsupported solvents

    Example:
        >>> get_alpb_solvent("CHCl3")
        'chcl3'
        >>> get_alpb_solvent("water")
        None
    """
    # Normalize to lowercase for case-insensitive matching
    solvent_lower = job_solvent.lower()

    # Return None for vacuum
    if solvent_lower == "vacuum":
        return None

    # Look up in ALPB mapping
    return ALPB_SOLVENT_MAP.get(solvent_lower)


def parse_crest_ensemble(ensemble_file: Path) -> list[CRESTConformer]:
    """
    Parse CREST multi-structure XYZ ensemble file.

    CREST outputs conformers in concatenated XYZ format:
    - Each structure starts with atom count line
    - Comment line contains energy as first token (in Hartree)
    - Coordinate lines follow
    - Structures are concatenated sequentially

    Args:
        ensemble_file: Path to CREST ensemble XYZ file

    Returns:
        List of CRESTConformer objects with sequential IDs, energies, and xyz_blocks

    Example:
        >>> conformers = parse_crest_ensemble(Path("crest_conformers.xyz"))
        >>> len(conformers)
        5
        >>> conformers[0].conformer_id
        'conf_001'
        >>> conformers[0].energy_hartree
        -12.345
    """
    conformers = []
    lines = ensemble_file.read_text().strip().split("\n")

    i = 0
    conf_number = 1

    while i < len(lines):
        # Read atom count
        if not lines[i].strip():
            i += 1
            continue

        num_atoms = int(lines[i].strip())

        # Read comment line (contains energy)
        i += 1
        if i >= len(lines):
            break

        comment_line = lines[i].strip()
        energy_hartree = float(comment_line.split()[0])

        # Read coordinate lines
        i += 1
        coord_lines = []
        for _ in range(num_atoms):
            if i >= len(lines):
                break
            coord_lines.append(lines[i])
            i += 1

        # Build xyz_block
        xyz_block = f"{num_atoms}\n{comment_line}\n"
        xyz_block += "\n".join(coord_lines) + "\n"

        # Create conformer with sequential ID
        conformer_id = f"conf_{conf_number:03d}"
        conformer = CRESTConformer(
            conformer_id=conformer_id,
            energy_hartree=energy_hartree,
            xyz_block=xyz_block,
        )
        conformers.append(conformer)

        conf_number += 1

    return conformers
