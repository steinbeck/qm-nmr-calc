"""NMReData SDF file generation for NMR prediction results.

NMReData is a standardized format for NMR chemical shift assignments embedded in
SDF files. This module generates NMReData-compliant SDF files from calculated
NMR shifts.

Specification: https://nmredata.org/wiki/NMReDATA_tag_format
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

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


def format_sdf_tag(name: str, value: str) -> str:
    """Format a single SDF tag in NMReData format.

    Args:
        name: Tag name (e.g., 'NMREDATA_VERSION').
        value: Tag value content.

    Returns:
        Formatted SDF tag block with proper spacing.

    Example:
        >>> format_sdf_tag("NMREDATA_VERSION", "1.1")
        '>  <NMREDATA_VERSION>\\n1.1\\n'
    """
    return f">  <{name}>\n{value}\n"


def generate_nmredata_sdf(
    smiles: str,
    geometry_xyz: str,
    h1_shifts: list[dict],
    c13_shifts: list[dict],
    solvent: str,
    temperature_k: float = 298.15,
    functional: str = "B3LYP",
    basis_set: str = "6-311+G(2d,p)",
    is_ensemble: bool = False,
    conformer_count: int | None = None,
) -> str:
    """Generate NMReData-formatted SDF file content.

    Args:
        smiles: Input SMILES string.
        geometry_xyz: XYZ coordinate file content (optimized geometry).
        h1_shifts: 1H chemical shifts with 1-based atom indices.
            List of dicts with keys: index, atom, shift.
        c13_shifts: 13C chemical shifts with 1-based atom indices.
            List of dicts with keys: index, atom, shift.
        solvent: NWChem COSMO solvent name (chcl3, dmso, vacuum).
        temperature_k: Temperature in Kelvin (for Boltzmann averaging).
        functional: DFT functional used (e.g., 'B3LYP').
        basis_set: Basis set used (e.g., '6-311+G(2d,p)').
        is_ensemble: Whether shifts are Boltzmann-averaged from ensemble.
        conformer_count: Number of conformers (for ensemble metadata).

    Returns:
        Complete NMReData SDF file content as string.

    Raises:
        ValueError: If SMILES is invalid or XYZ atom count mismatches.

    Example:
        >>> sdf = generate_nmredata_sdf(
        ...     smiles="CCO",
        ...     geometry_xyz=xyz_content,
        ...     h1_shifts=[{"index": 4, "atom": "H", "shift": 1.18}],
        ...     c13_shifts=[{"index": 1, "atom": "C", "shift": 18.2}],
        ...     solvent="chcl3"
        ... )
        >>> "NMREDATA_VERSION" in sdf
        True
    """
    # Parse XYZ to extract coordinates
    xyz_lines = geometry_xyz.strip().split("\n")
    if len(xyz_lines) < 3:
        raise ValueError("Invalid XYZ format: need at least 3 lines")

    # First line is atom count, second is comment, rest are coords
    try:
        atom_count = int(xyz_lines[0].strip())
    except ValueError:
        raise ValueError(f"Invalid XYZ: first line must be atom count, got '{xyz_lines[0]}'")

    coords = []
    for line in xyz_lines[2:]:
        parts = line.split()
        if len(parts) >= 4:
            # XYZ format: element x y z
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

    if len(coords) != atom_count:
        raise ValueError(
            f"XYZ atom count mismatch: header says {atom_count}, found {len(coords)} coordinates"
        )

    # Create RDKit Mol from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: '{smiles}'")

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Verify atom count matches
    if mol.GetNumAtoms() != atom_count:
        raise ValueError(
            f"Atom count mismatch: RDKit mol has {mol.GetNumAtoms()} atoms, "
            f"XYZ has {atom_count} atoms"
        )

    # Create conformer and set 3D coordinates
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        if i < mol.GetNumAtoms():
            conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf, assignId=True)

    # Generate base MOL block
    mol_block = Chem.MolToMolBlock(mol)

    # Build NMReData tags
    tags = []

    # Required tags
    tags.append(format_sdf_tag("NMREDATA_VERSION", NMREDATA_VERSION))
    tags.append(format_sdf_tag("NMREDATA_LEVEL", NMREDATA_LEVEL))
    tags.append(format_sdf_tag("NMREDATA_SOLVENT", map_solvent_to_nmredata(solvent)))
    tags.append(format_sdf_tag("NMREDATA_TEMPERATURE", f"{temperature_k:.2f}"))

    # ASSIGNMENT tag - core chemical shift data
    assignment_content = format_assignment_tag(h1_shifts, c13_shifts)
    if assignment_content:  # Only add if non-empty
        tags.append(format_sdf_tag("NMREDATA_ASSIGNMENT", assignment_content))

    # Optional but recommended tags
    formula = rdMolDescriptors.CalcMolFormula(mol)
    tags.append(format_sdf_tag("NMREDATA_FORMULA", formula))
    tags.append(format_sdf_tag("NMREDATA_SMILES", smiles))

    # Provenance tag - document prediction method
    provenance_lines = [
        "Predicted by qm-nmr-calc",
        f"Method: {functional}/{basis_set}",
        "Scaling: DELTA50",
    ]
    if is_ensemble and conformer_count:
        provenance_lines.append(
            f"Boltzmann-averaged from {conformer_count} conformers at {temperature_k} K"
        )
    provenance = "\n".join(provenance_lines)
    tags.append(format_sdf_tag("NMREDATA_ID", provenance))

    # Assemble complete SDF: MOL block + tags + terminator
    sdf_content = mol_block
    for tag in tags:
        sdf_content += "\n" + tag

    # SDF terminator
    sdf_content += "\n$$$$\n"

    return sdf_content
