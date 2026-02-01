"""NMReData SDF file generation for NMR prediction results.

NMReData is a standardized format for NMR chemical shift assignments embedded in
SDF files. This module generates NMReData-compliant SDF files from calculated
NMR shifts.

Specification: https://nmredata.org/wiki/NMReDATA_tag_format
"""

# NMReData format constants
NMREDATA_VERSION = "1.1"  # Stable NMReData format version
NMREDATA_LEVEL = "0"  # Level 0: predicted data with unambiguous assignments
NMREDATA_SEP = ", "  # Standard separator: comma + space


def map_solvent_to_nmredata(solvent: str) -> str:
    """Map NWChem COSMO solvent names to NMReData solvent names.

    Args:
        solvent: NWChem COSMO solvent name (lowercase, from our system).

    Returns:
        NMReData solvent name (deuterated form for NMR solvents).

    Raises:
        ValueError: If solvent is not recognized.

    Examples:
        >>> map_solvent_to_nmredata("chcl3")
        'CDCl3'
        >>> map_solvent_to_nmredata("dmso")
        '(CD3)2SO'
        >>> map_solvent_to_nmredata("vacuum")
        'vacuum'
    """
    # Mapping from NWChem names to NMReData convention (deuterated forms)
    solvent_map = {
        "chcl3": "CDCl3",
        "dmso": "(CD3)2SO",  # NMReData convention for deuterated DMSO
        "vacuum": "vacuum",  # Gas phase - no solvent
    }

    normalized = solvent.lower().strip()
    if normalized not in solvent_map:
        raise ValueError(
            f"Unknown solvent '{solvent}'. Supported: {list(solvent_map.keys())}"
        )

    return solvent_map[normalized]


def format_atom_label(atom: str, index: int) -> str:
    """Format atom label for NMReData ASSIGNMENT tag.

    Args:
        atom: Element symbol ('H' or 'C').
        index: 1-based atom index.

    Returns:
        Lowercase label in format: element + index (e.g., 'h5', 'c3').

    Examples:
        >>> format_atom_label("H", 5)
        'h5'
        >>> format_atom_label("C", 3)
        'c3'
    """
    return f"{atom.lower()}{index}"
