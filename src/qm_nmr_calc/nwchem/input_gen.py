"""NWChem input file generation for geometry optimization and NMR shielding.

Input format follows ISiCLE conventions for compatibility and reliability.
"""

# Supported solvents for COSMO (NWChem recognizes these by name)
# "vacuum" is special - no COSMO block will be generated for gas-phase calculations
SUPPORTED_SOLVENTS = {"chcl3", "dmso", "water", "acetone", "methanol", "benzene", "vacuum"}


def _validate_solvent(solvent: str) -> str:
    """Validate solvent and return normalized name.

    Args:
        solvent: Solvent name (case-insensitive)

    Returns:
        Normalized solvent name for NWChem COSMO

    Raises:
        ValueError: If solvent is not supported
    """
    solvent_lower = solvent.lower()
    if solvent_lower not in SUPPORTED_SOLVENTS:
        valid_solvents = ", ".join(sorted(SUPPORTED_SOLVENTS))
        raise ValueError(
            f"Unsupported solvent: '{solvent}'. "
            f"Valid options: {valid_solvents}"
        )
    return solvent_lower


def generate_optimization_input(
    geometry_xyz: str,
    functional: str,
    basis_set: str,
    solvent: str,
    max_iter: int = 150,
) -> str:
    """Generate NWChem input file for geometry optimization.

    Uses ISiCLE-compatible format:
    - Memory allocation (1600 MB global, 100 MB heap, 600 MB stack)
    - noautoz noautosym geometry options
    - Driver block for optimization control
    - COSMO solvation block (omitted for vacuum/gas-phase calculations)

    Args:
        geometry_xyz: XYZ-format geometry (atom lines only, no header)
        functional: DFT functional (e.g., 'b3lyp')
        basis_set: Basis set name (e.g., '6-31G*')
        solvent: Solvent name for COSMO (chcl3, dmso, etc.) or "vacuum" for gas-phase
        max_iter: Maximum optimization iterations

    Returns:
        Complete NWChem input file as string

    Raises:
        ValueError: If solvent is not supported
    """
    solvent_name = _validate_solvent(solvent)

    # Build COSMO block only for non-vacuum calculations
    if solvent_name == "vacuum":
        cosmo_block = ""
    else:
        cosmo_block = f"""
cosmo
  do_gasphase False
  solvent {solvent_name}
end
"""

    return f"""start molecule
title "Geometry Optimization"

memory global 1600 mb heap 100 mb stack 600 mb

geometry units angstrom noautoz noautosym
{geometry_xyz}
end

basis spherical
  * library {basis_set}
end

dft
  xc {functional}
end

driver
  maxiter {max_iter}
  xyz molecule_geom
end
{cosmo_block}
task dft optimize
"""


def generate_shielding_input(
    geometry_xyz: str,
    functional: str,
    basis_set: str,
    solvent: str,
) -> str:
    """Generate NWChem input file for NMR shielding calculation.

    Uses ISiCLE-compatible format. COSMO solvation is omitted for vacuum/gas-phase.

    The DFT block includes 'direct' to force on-the-fly integral evaluation,
    which is required for stable CPHF (coupled-perturbed Hartree-Fock) convergence
    during NMR property calculations with COSMO solvation. Without 'direct',
    NWChem may fail with 'cphf_solve2: SCF residual greater than 1d-2'.

    Args:
        geometry_xyz: XYZ-format geometry (atom lines only, no header)
        functional: DFT functional (e.g., 'b3lyp')
        basis_set: Basis set name (e.g., '6-311+G(2d,p)')
        solvent: Solvent name for COSMO (chcl3, dmso, etc.) or "vacuum" for gas-phase

    Returns:
        Complete NWChem input file as string

    Raises:
        ValueError: If solvent is not supported
    """
    solvent_name = _validate_solvent(solvent)

    # Build COSMO block only for non-vacuum calculations
    if solvent_name == "vacuum":
        cosmo_block = ""
    else:
        cosmo_block = f"""
cosmo
  do_gasphase False
  solvent {solvent_name}
end
"""

    return f"""start molecule
title "NMR Shielding Calculation"

memory global 1600 mb heap 100 mb stack 600 mb

geometry units angstrom noautoz noautosym
{geometry_xyz}
end

basis spherical
  * library {basis_set}
end

dft
  xc {functional}
  direct
end
{cosmo_block}
property
  shielding
end

task dft property
"""
