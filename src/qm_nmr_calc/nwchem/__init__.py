"""NWChem input/output handling for NMR calculations.

This module provides direct NWChem integration for geometry optimization
and NMR shielding calculations with COSMO solvation, replacing the previous
ISiCLE dependency.
"""

from qm_nmr_calc.nwchem.geometry import (
    load_geometry_file,
    mol_to_xyz_block,
    smiles_to_xyz,
    validate_geometry,
)
from qm_nmr_calc.nwchem.input_gen import (
    generate_optimization_input,
    generate_shielding_input,
)

__all__ = [
    # Geometry handling
    "smiles_to_xyz",
    "load_geometry_file",
    "mol_to_xyz_block",
    "validate_geometry",
    # Input generation
    "generate_optimization_input",
    "generate_shielding_input",
]
