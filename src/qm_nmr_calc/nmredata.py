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


def format_assignment_tag(h1_shifts: list, c13_shifts: list) -> str:
    """Format NMREDATA_ASSIGNMENT tag content from shift data.

    Assumes input indices are 1-based (NWChem convention). Does NOT perform
    any index conversion - uses indices as-is from the input data.

    Args:
        h1_shifts: List of 1H shift dicts with keys: index (1-based), atom, shift.
        c13_shifts: List of 13C shift dicts with keys: index (1-based), atom, shift.

    Returns:
        Multiline ASSIGNMENT tag content with format: label, shift, index.
        Returns empty string if both lists are empty.

    Format per line:
        {label}, {shift:.4f}, {atom_index}

    Example:
        >>> h1 = [{"index": 1, "atom": "H", "shift": 7.2453}]
        >>> c13 = [{"index": 2, "atom": "C", "shift": 128.45}]
        >>> format_assignment_tag(h1, c13)
        'h1, 7.2453, 1\\nc2, 128.4500, 2'
    """
    lines = []

    # Process 1H shifts
    for shift_data in h1_shifts:
        label = format_atom_label(shift_data["atom"], shift_data["index"])
        shift_value = shift_data["shift"]
        atom_index = shift_data["index"]
        # Format: label, shift (4 decimals), atom_index
        line = f"{label}{NMREDATA_SEP}{shift_value:.4f}{NMREDATA_SEP}{atom_index}"
        lines.append(line)

    # Process 13C shifts
    for shift_data in c13_shifts:
        label = format_atom_label(shift_data["atom"], shift_data["index"])
        shift_value = shift_data["shift"]
        atom_index = shift_data["index"]
        # Format: label, shift (4 decimals), atom_index
        line = f"{label}{NMREDATA_SEP}{shift_value:.4f}{NMREDATA_SEP}{atom_index}"
        lines.append(line)

    return "\n".join(lines)
