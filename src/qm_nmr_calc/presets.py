"""Calculation preset configurations for NMR calculations."""

import os
from enum import Enum
from typing import TypedDict


def get_default_processes() -> int:
    """Get default number of MPI processes from environment or CPU count.

    Uses NWCHEM_NPROC env var if set, otherwise detects available CPUs
    and caps at 40 (diminishing returns beyond that for typical molecules).
    """
    env_val = os.environ.get("NWCHEM_NPROC")
    if env_val:
        try:
            return int(env_val)
        except ValueError:
            pass

    # Detect available CPUs, cap at 40
    try:
        cpu_count = os.cpu_count() or 4
        return min(cpu_count, 40)
    except Exception:
        return 4


class PresetName(str, Enum):
    """Available calculation presets."""

    DRAFT = "draft"
    PRODUCTION = "production"


class CalculationPreset(TypedDict):
    """Configuration for a calculation preset.

    Attributes:
        name: Human-readable preset name
        description: Description of the preset's purpose
        functional: DFT functional (e.g., 'b3lyp')
        basis_set: Basis set for geometry optimization
        nmr_basis_set: Basis set for NMR shielding calculation
        processes: Number of parallel processes
        max_iter: Maximum geometry optimization iterations
    """

    name: str
    description: str
    functional: str
    basis_set: str
    nmr_basis_set: str
    processes: int
    max_iter: int


# Compute default processes once at module load
_DEFAULT_PROCESSES = get_default_processes()

PRESETS: dict[PresetName, CalculationPreset] = {
    PresetName.DRAFT: {
        "name": "draft",
        "description": "Fast calculations for quick checks",
        "functional": "b3lyp",
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-31G*",
        "processes": _DEFAULT_PROCESSES,
        "max_iter": 100,
    },
    PresetName.PRODUCTION: {
        "name": "production",
        "description": "Balanced accuracy and compute time (default)",
        "functional": "b3lyp",
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-311+G(2d,p)",
        "processes": _DEFAULT_PROCESSES,
        "max_iter": 300,
    },
}

DEFAULT_PRESET = PresetName.PRODUCTION
