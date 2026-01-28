"""
CREST conformer generation utilities.

This module provides utilities for high-accuracy conformer generation using
CREST/xTB when available:
- Binary detection for CREST and xTB
- ALPB solvent model mapping
- Multi-structure XYZ file parsing
- CREST subprocess execution
- Ensemble generation pipeline
"""

import os
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple

from qm_nmr_calc.conformers.boltzmann import HARTREE_TO_KCAL
from qm_nmr_calc.conformers.filters import filter_by_energy_window
from qm_nmr_calc.models import ConformerData, ConformerEnsemble
from qm_nmr_calc.nwchem.geometry import smiles_to_xyz
from qm_nmr_calc.storage import (
    create_conformer_directories,
    get_conformer_output_dir,
    get_job_dir,
)


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


def run_crest(
    input_xyz: Path,
    solvent: str,
    charge: int = 0,
    ewin_kcal: float = 6.0,
    num_threads: int | None = None,
    timeout_seconds: int = 7200,
) -> Path:
    """
    Run CREST conformational search subprocess.

    Executes CREST with GFN2-xTB and ALPB implicit solvent. Includes timeout
    handling for macrocycles that may hang.

    Args:
        input_xyz: Path to input XYZ file
        solvent: ALPB solvent name (e.g., "chcl3", "dmso")
        charge: Molecular charge (default: 0)
        ewin_kcal: Energy window in kcal/mol for conformer filtering (default: 6.0)
        num_threads: Number of threads (default: os.cpu_count() or 4)
        timeout_seconds: Subprocess timeout in seconds (default: 7200 = 2 hours)

    Returns:
        Path to crest_conformers.xyz output file

    Raises:
        RuntimeError: If CREST times out, exits with error, or output file missing

    Example:
        >>> output = run_crest(Path("input.xyz"), "chcl3", charge=0)
        >>> output.exists()
        True
    """
    # Default num_threads
    if num_threads is None:
        num_threads = os.cpu_count() or 4

    # Build CREST command
    cmd = [
        "crest",
        str(input_xyz),
        "--gfn2",
        "--alpb",
        solvent,
        "--chrg",
        str(charge),
        "--ewin",
        str(ewin_kcal),
        "-T",
        str(num_threads),
    ]

    # Set up environment variables
    env = os.environ.copy()
    env["OMP_STACKSIZE"] = "2G"
    env["GFORTRAN_UNBUFFERED_ALL"] = "1"

    # Run CREST subprocess
    try:
        subprocess.run(
            cmd,
            cwd=input_xyz.parent,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            env=env,
            check=True,
        )
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"CREST timeout after {timeout_seconds} seconds. "
            "This can happen with large/complex molecules (especially macrocycles). "
            "Consider using RDKit conformer generation mode instead."
        ) from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"CREST failed with exit code {e.returncode}. "
            f"Error output: {e.stderr}"
        ) from e

    # Verify output file exists
    output_file = input_xyz.parent / "crest_conformers.xyz"
    if not output_file.exists():
        raise RuntimeError(
            f"CREST completed but output file not found: {output_file}"
        )

    return output_file


def generate_crest_ensemble(
    smiles: str,
    job_id: str,
    solvent: str,
    charge: int = 0,
    energy_window_kcal: float = 6.0,
    timeout_seconds: int = 7200,
) -> ConformerEnsemble:
    """
    Generate conformer ensemble using CREST.

    Pipeline:
    1. Generate initial 3D geometry from SMILES using RDKit
    2. Run CREST conformational search
    3. Parse CREST output
    4. Convert energies from Hartree to relative kcal/mol
    5. Apply energy window filter
    6. Write individual XYZ files
    7. Build ConformerEnsemble

    Args:
        smiles: SMILES string
        job_id: Job identifier
        solvent: ALPB solvent name (must be pre-validated)
        charge: Molecular charge (default: 0)
        energy_window_kcal: Energy window for filtering (default: 6.0 kcal/mol)
        timeout_seconds: CREST timeout (default: 7200 = 2 hours)

    Returns:
        ConformerEnsemble with method="crest", relative energies in kcal/mol

    Example:
        >>> ensemble = generate_crest_ensemble("C", "job123", "chcl3")
        >>> ensemble.method
        'crest'
        >>> len(ensemble.conformers) > 0
        True
    """
    # Get job directory and create input XYZ
    job_dir = get_job_dir(job_id)
    output_dir = job_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate initial geometry from SMILES
    mol, xyz_block = smiles_to_xyz(smiles)
    input_xyz = output_dir / "crest_input.xyz"
    input_xyz.write_text(xyz_block)

    # Run CREST
    crest_output = run_crest(
        input_xyz,
        solvent,
        charge=charge,
        ewin_kcal=energy_window_kcal,
        timeout_seconds=timeout_seconds,
    )

    # Parse CREST output
    crest_conformers = parse_crest_ensemble(crest_output)
    total_generated = len(crest_conformers)

    # Convert energies from Hartree to relative kcal/mol
    hartree_energies = [conf.energy_hartree for conf in crest_conformers]
    min_energy_hartree = min(hartree_energies)

    # Calculate relative energies in kcal/mol
    relative_energies_kcal = [
        (energy - min_energy_hartree) * HARTREE_TO_KCAL
        for energy in hartree_energies
    ]

    # Apply energy window filter
    # filter_by_energy_window expects conf_ids as list[int], but we'll use indices
    conf_indices = list(range(len(crest_conformers)))
    filtered_indices, filtered_energies = filter_by_energy_window(
        conf_indices, relative_energies_kcal, window_kcal=energy_window_kcal
    )

    # Filter conformers
    filtered_conformers = [crest_conformers[i] for i in filtered_indices]
    total_after_pre_filter = len(filtered_conformers)

    # Create per-conformer directories
    conformer_ids = [conf.conformer_id for conf in filtered_conformers]
    create_conformer_directories(job_id, conformer_ids)

    # Write individual XYZ files and build ConformerData list
    conformer_data_list = []
    for conf, energy_kcal in zip(filtered_conformers, filtered_energies):
        # Get output directory for this conformer
        conf_output_dir = get_conformer_output_dir(job_id, conf.conformer_id)

        # Ensure directory exists (create_conformer_directories should have done this)
        conf_output_dir.mkdir(parents=True, exist_ok=True)

        # Write XYZ file
        xyz_file = conf_output_dir / "geometry.xyz"
        xyz_file.write_text(conf.xyz_block)

        # Build ConformerData
        conformer_data = ConformerData(
            conformer_id=conf.conformer_id,
            energy=energy_kcal,
            energy_unit="kcal_mol",
            geometry_file=f"output/conformers/{conf.conformer_id}/geometry.xyz",
            status="pending",
        )
        conformer_data_list.append(conformer_data)

    # Build ConformerEnsemble
    ensemble = ConformerEnsemble(
        method="crest",
        conformers=conformer_data_list,
        pre_dft_energy_window_kcal=energy_window_kcal,
        total_generated=total_generated,
        total_after_pre_filter=total_after_pre_filter,
    )

    return ensemble
