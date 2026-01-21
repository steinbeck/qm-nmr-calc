"""NWChem input/output handling for NMR calculations.

This module provides direct NWChem integration for geometry optimization
and NMR shielding calculations with COSMO solvation, replacing the previous
ISiCLE dependency.
"""

from qm_nmr_calc.nwchem.input_gen import (
    generate_optimization_input,
    generate_shielding_input,
)

__all__ = [
    "generate_optimization_input",
    "generate_shielding_input",
]
