"""NWChem COSMO solvent validation for NMR calculations."""

# Supported NMR solvents for COSMO solvation model.
# Keys are NWChem-compatible names (lowercase), values are descriptions.
# Limited to solvents with known dielectric constants in input_gen.py.
SUPPORTED_SOLVENTS: dict[str, str] = {
    "chcl3": "Chloroform (CDCl3)",
    "dmso": "Dimethylsulfoxide (DMSO-d6)",
    "vacuum": "Gas phase (no solvent)",
    "methanol": "Methanol (Methanol-d4)",
    "water": "Water (D2O)",
    "acetone": "Acetone (Acetone-d6)",
    "benzene": "Benzene (Benzene-d6)",
    "pyridine": "Pyridine (Pyridine-d5)",
    "thf": "Tetrahydrofuran (THF-d8)",
    "toluene": "Toluene (Toluene-d8)",
    "dcm": "Dichloromethane (CD2Cl2)",
    "acetonitrile": "Acetonitrile (CD3CN)",
    "dmf": "N,N-Dimethylformamide (DMF-d7)",
}


def validate_solvent(solvent: str) -> str | None:
    """Validate solvent name against supported NWChem COSMO solvents.

    Args:
        solvent: User-provided solvent name.

    Returns:
        Normalized solvent name if valid, None if invalid.
    """
    normalized = solvent.lower().strip()
    if normalized in SUPPORTED_SOLVENTS:
        return normalized
    return None


def get_supported_solvents() -> list[str]:
    """Get list of supported solvent names for API documentation.

    Returns:
        Sorted list of supported solvent names.
    """
    return sorted(SUPPORTED_SOLVENTS.keys())


def get_solvent_display_name(solvent: str) -> str:
    """Get display name for a solvent (deuterated NMR form).

    Args:
        solvent: NWChem solvent name (e.g., 'chcl3').

    Returns:
        Deuterated display name (e.g., 'CDCl3'), or original if not found.
    """
    desc = SUPPORTED_SOLVENTS.get(solvent.lower(), "")
    # Extract text in parentheses: "Chloroform (CDCl3)" -> "CDCl3"
    if "(" in desc and ")" in desc:
        start = desc.index("(") + 1
        end = desc.index(")")
        return desc[start:end]
    return solvent
