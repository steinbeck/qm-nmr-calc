"""NWChem input file generation for geometry optimization and NMR shielding."""

# COSMO dielectric constants for supported solvents
# Sources: Standard literature values
COSMO_DIELECTRIC: dict[str, float] = {
    "chcl3": 4.8,
    "dmso": 46.0,
}


def _validate_solvent(solvent: str) -> float:
    """Validate solvent and return its dielectric constant.

    Args:
        solvent: Solvent name (case-insensitive)

    Returns:
        Dielectric constant for the solvent

    Raises:
        ValueError: If solvent is not supported
    """
    solvent_lower = solvent.lower()
    if solvent_lower not in COSMO_DIELECTRIC:
        valid_solvents = ", ".join(sorted(COSMO_DIELECTRIC.keys()))
        raise ValueError(
            f"Unsupported solvent: '{solvent}'. "
            f"Valid options: {valid_solvents}"
        )
    return COSMO_DIELECTRIC[solvent_lower]


def generate_optimization_input(
    geometry_xyz: str,
    functional: str,
    basis_set: str,
    solvent: str,
    max_iter: int = 150,
) -> str:
    """Generate NWChem input file for geometry optimization with COSMO.

    Args:
        geometry_xyz: XYZ-format geometry (atom lines only, no header)
        functional: DFT functional (e.g., 'b3lyp')
        basis_set: Basis set name (e.g., '6-31G*')
        solvent: Solvent name for COSMO (chcl3 or dmso)
        max_iter: Maximum optimization iterations

    Returns:
        Complete NWChem input file as string

    Raises:
        ValueError: If solvent is not supported
    """
    dielec = _validate_solvent(solvent)

    return f"""start molecule
title "Geometry Optimization"

geometry units angstrom noautosym
{geometry_xyz}
end

basis spherical
  * library "{basis_set}"
end

dft
  xc {functional}
  iterations {max_iter}
  convergence energy 1e-7 density 1e-5 gradient 5e-4
  direct
end

cosmo
  dielec {dielec}
end

task dft optimize
"""


def generate_shielding_input(
    geometry_xyz: str,
    functional: str,
    basis_set: str,
    solvent: str,
) -> str:
    """Generate NWChem input file for NMR shielding calculation with COSMO.

    Args:
        geometry_xyz: XYZ-format geometry (atom lines only, no header)
        functional: DFT functional (e.g., 'b3lyp')
        basis_set: Basis set name (e.g., '6-311+G(2d,p)')
        solvent: Solvent name for COSMO (chcl3 or dmso)

    Returns:
        Complete NWChem input file as string

    Raises:
        ValueError: If solvent is not supported
    """
    dielec = _validate_solvent(solvent)

    return f"""start molecule
title "NMR Shielding Calculation"

geometry units angstrom noautosym
{geometry_xyz}
end

basis spherical
  * library "{basis_set}"
end

dft
  xc {functional}
  direct
end

cosmo
  dielec {dielec}
end

property
  shielding
end

task dft property
"""
