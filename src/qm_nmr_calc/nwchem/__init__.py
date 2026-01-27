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
from qm_nmr_calc.nwchem.output_parser import (
    extract_dft_energy,
    extract_optimized_geometry,
    parse_shielding_output,
)
from qm_nmr_calc.nwchem.runner import (
    apply_post_dft_filter,
    get_nwchem_version,
    run_calculation,
    run_conformer_dft_optimization,
    run_conformer_nmr_calculations,
    run_ensemble_dft_and_nmr,
    run_nwchem,
    validate_nwchem,
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
    # Output parsing
    "extract_dft_energy",
    "extract_optimized_geometry",
    "parse_shielding_output",
    # Runner
    "run_nwchem",
    "run_calculation",
    "run_conformer_dft_optimization",
    "run_conformer_nmr_calculations",
    "run_ensemble_dft_and_nmr",
    "apply_post_dft_filter",
    "validate_nwchem",
    "get_nwchem_version",
]
