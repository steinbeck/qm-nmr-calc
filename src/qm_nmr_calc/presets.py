"""Calculation preset configurations for NMR calculations."""

from enum import Enum
from typing import TypedDict


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


PRESETS: dict[PresetName, CalculationPreset] = {
    PresetName.DRAFT: {
        "name": "draft",
        "description": "Fast calculations for quick checks",
        "functional": "b3lyp",
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-31G*",
        "processes": 40,
        "max_iter": 100,
    },
    PresetName.PRODUCTION: {
        "name": "production",
        "description": "Balanced accuracy and compute time (default)",
        "functional": "b3lyp",
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-311+G(2d,p)",
        "processes": 40,
        "max_iter": 300,
    },
}

DEFAULT_PRESET = PresetName.PRODUCTION
