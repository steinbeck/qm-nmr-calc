"""xTB-based conformer energy ranking for pre-DFT selection.

GFN2-xTB provides fast semi-empirical energies (100-1000x faster than DFT)
with better ranking correlation than MMFF force fields.

Requires xTB binary in PATH. Install via: conda install -c conda-forge xtb
"""

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem import rdmolfiles


# Conversion factor
HARTREE_TO_KCAL = 627.509474


def detect_xtb_available() -> bool:
    """Check if xTB binary is available in PATH.

    Returns:
        True if xTB is found and executable, False otherwise.
    """
    return shutil.which("xtb") is not None


def get_xtb_version() -> Optional[str]:
    """Get xTB version string if available.

    Returns:
        Version string (e.g., "6.6.1") or None if xTB not available.
    """
    if not detect_xtb_available():
        return None

    try:
        result = subprocess.run(
            ["xtb", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        # Parse version from output
        for line in result.stdout.split("\n"):
            if "version" in line.lower():
                parts = line.split()
                for i, part in enumerate(parts):
                    if part.lower() == "version" and i + 1 < len(parts):
                        return parts[i + 1]
        return "unknown"
    except Exception:
        return None


def calculate_xtb_energy(
    xyz_content: str,
    charge: int = 0,
    multiplicity: int = 1,
    solvent: Optional[str] = None,
    timeout_seconds: int = 60,
) -> float:
    """Calculate GFN2-xTB single-point energy for a conformer.

    Args:
        xyz_content: XYZ file content as string
        charge: Molecular charge (default 0)
        multiplicity: Spin multiplicity (default 1 for singlet)
        solvent: ALPB solvent name (e.g., "chcl3", "dmso") or None for gas phase
        timeout_seconds: Maximum time for calculation (default 60s)

    Returns:
        Total energy in Hartree

    Raises:
        RuntimeError: If xTB not available or calculation fails
        TimeoutError: If calculation exceeds timeout
    """
    if not detect_xtb_available():
        raise RuntimeError("xTB binary not found in PATH. Install via: conda install -c conda-forge xtb")

    with tempfile.TemporaryDirectory() as tmpdir:
        xyz_path = Path(tmpdir) / "input.xyz"
        xyz_path.write_text(xyz_content)

        # Build command
        cmd = [
            "xtb",
            str(xyz_path),
            "--gfn2",           # Use GFN2-xTB method
            "--sp",             # Single point (no optimization)
            "--chrg", str(charge),
            "--uhf", str(multiplicity - 1),  # xTB uses number of unpaired electrons
        ]

        # Add solvation if specified
        if solvent:
            alpb_solvent = _map_solvent_to_alpb(solvent)
            if alpb_solvent:
                cmd.extend(["--alpb", alpb_solvent])

        try:
            result = subprocess.run(
                cmd,
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
            )
        except subprocess.TimeoutExpired:
            raise TimeoutError(f"xTB calculation timed out after {timeout_seconds}s")

        if result.returncode != 0:
            raise RuntimeError(f"xTB failed with code {result.returncode}: {result.stderr[:500]}")

        # Parse energy from output
        return _parse_xtb_energy(result.stdout)


def _map_solvent_to_alpb(solvent: str) -> Optional[str]:
    """Map common solvent names to xTB ALPB solvent names.

    Args:
        solvent: Solvent name (various formats accepted)

    Returns:
        ALPB solvent name or None if not supported
    """
    solvent_map = {
        # Chloroform variants
        "chcl3": "chcl3",
        "chloroform": "chcl3",
        "cdcl3": "chcl3",
        # DMSO variants
        "dmso": "dmso",
        "dmso-d6": "dmso",
        # Water
        "water": "water",
        "h2o": "water",
        "d2o": "water",
        # Methanol
        "methanol": "methanol",
        "meoh": "methanol",
        "cd3od": "methanol",
        # Acetone
        "acetone": "acetone",
        # Acetonitrile
        "acetonitrile": "acetonitrile",
        "mecn": "acetonitrile",
        # THF
        "thf": "thf",
        # Benzene
        "benzene": "benzene",
        "c6d6": "benzene",
        # Toluene
        "toluene": "toluene",
        # DCM
        "dcm": "ch2cl2",
        "ch2cl2": "ch2cl2",
        "dichloromethane": "ch2cl2",
    }
    return solvent_map.get(solvent.lower())


def _parse_xtb_energy(output: str) -> float:
    """Parse total energy from xTB output.

    Args:
        output: xTB stdout content

    Returns:
        Energy in Hartree

    Raises:
        RuntimeError: If energy cannot be parsed
    """
    # Look for "TOTAL ENERGY" line
    # Format: "          | TOTAL ENERGY              -XX.XXXXXX Eh   |"
    for line in output.split("\n"):
        if "TOTAL ENERGY" in line and "Eh" in line:
            parts = line.split()
            for i, part in enumerate(parts):
                if part == "ENERGY" and i + 1 < len(parts):
                    try:
                        return float(parts[i + 1])
                    except ValueError:
                        continue

    raise RuntimeError("Could not parse energy from xTB output")


def rank_conformers_by_xtb(
    mol: Chem.Mol,
    conf_ids: list[int],
    charge: int = 0,
    solvent: Optional[str] = None,
    timeout_per_conf: int = 60,
) -> dict[int, float]:
    """Calculate xTB energies for multiple conformers.

    Args:
        mol: RDKit Mol with conformers
        conf_ids: List of conformer IDs to rank
        charge: Molecular charge
        solvent: Solvent name or None for gas phase
        timeout_per_conf: Timeout in seconds per conformer

    Returns:
        Dict mapping conformer ID to relative energy in kcal/mol.
        Energies are relative to the minimum (lowest = 0.0).

    Raises:
        RuntimeError: If all conformers fail or xTB not available
    """
    if not detect_xtb_available():
        raise RuntimeError("xTB not available")

    energies_hartree = {}
    failures = []

    for conf_id in conf_ids:
        # Generate XYZ content for this conformer
        xyz_content = rdmolfiles.MolToXYZBlock(mol, confId=conf_id)

        try:
            energy = calculate_xtb_energy(
                xyz_content,
                charge=charge,
                solvent=solvent,
                timeout_seconds=timeout_per_conf,
            )
            energies_hartree[conf_id] = energy
        except Exception as e:
            # Log warning but continue with other conformers
            failures.append((conf_id, str(e)))
            continue

    if not energies_hartree:
        failure_msgs = "; ".join(f"{cid}: {msg}" for cid, msg in failures[:3])
        raise RuntimeError(f"All xTB calculations failed. First failures: {failure_msgs}")

    # Convert to kcal/mol relative to minimum
    min_energy = min(energies_hartree.values())

    return {
        conf_id: (energy - min_energy) * HARTREE_TO_KCAL
        for conf_id, energy in energies_hartree.items()
    }
