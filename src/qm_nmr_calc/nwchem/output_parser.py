"""Parse NWChem output files for geometry and NMR shielding data.

This module extracts optimized geometries and isotropic shielding tensors
from NWChem DFT calculation output files.

Expected NWChem Output Formats (NWChem 7.x):

Geometry Optimization:
    The optimized geometry appears after "Output coordinates in angstroms"
    with columns: No. Tag Charge X Y Z
    Example:
        Output coordinates in angstroms (scale by  1.889725989 to convert to a.u.)

          No.       Tag         Charge          X              Y              Z
         ---- ---------------- ---------- -------------- -------------- --------------
            1 C                    6.0000     0.00000000     0.00000000     0.00000000
            2 H                    1.0000     1.09000000     0.00000000     0.00000000

NMR Shielding:
    Shielding tensors appear under "GIAO Chemical Shielding Tensors"
    Each atom block shows the isotropic shielding value.
    Example:
        GIAO Chemical Shielding Tensors
        --------------------------------
        Atom:    1  C
           isotropic  =     183.4567
           ...

Note: These patterns are based on NWChem documentation and may need
adjustment for different NWChem versions. Test with actual output files.
"""

import re


def extract_optimized_geometry(output_text: str) -> str:
    """Extract the final optimized geometry from NWChem output.

    Searches for the "Output coordinates in angstroms" section and extracts
    atom coordinates in XYZ format (without header lines).

    Args:
        output_text: Full NWChem output file content as string.

    Returns:
        XYZ-format coordinate block as string, with lines like:
        "C     0.000000     0.000000     0.000000"

    Raises:
        RuntimeError: If no geometry section is found in the output.

    Example:
        >>> output = Path("optimize.out").read_text()
        >>> xyz = extract_optimized_geometry(output)
        >>> # Returns: "C     0.000000     0.000000     0.000000\\nH     1.090000..."
    """
    # Pattern to find the coordinate section header
    # NWChem prints "Output coordinates in angstroms" followed by a table
    coord_pattern = re.compile(
        r"Output coordinates in angstroms.*?\n"  # Header line
        r".*?No\.\s+Tag\s+Charge\s+X\s+Y\s+Z.*?\n"  # Column headers
        r"\s*-+.*?\n"  # Separator line (dashes)
        r"(.*?)(?:\n\s*\n|\Z)",  # Atom lines until blank line or end
        re.DOTALL | re.IGNORECASE,
    )

    match = coord_pattern.search(output_text)

    if not match:
        # Try alternative pattern for different NWChem versions
        # Some versions use "Geometry after optimization" or similar
        alt_pattern = re.compile(
            r"(?:Geometry|Output)\s+(?:after|coordinates).*?angstrom.*?\n"
            r".*?(?:atom|no\.)\s+.*?[xyz].*?\n"  # Header with coordinates
            r"\s*-+.*?\n"  # Separator
            r"(.*?)(?:\n\s*\n|\Z)",
            re.DOTALL | re.IGNORECASE,
        )
        match = alt_pattern.search(output_text)

    if not match:
        raise RuntimeError(
            "Could not find optimized geometry coordinates in NWChem output. "
            "Expected section starting with 'Output coordinates in angstroms'."
        )

    # Parse atom lines from the matched section
    coord_section = match.group(1)
    xyz_lines = []

    # Pattern for atom line: index, tag (element), charge, x, y, z
    # Example: "    1 C                    6.0000     0.00000000     0.00000000     0.00000000"
    atom_line_pattern = re.compile(
        r"^\s*\d+\s+"  # Index number
        r"([A-Za-z][a-z]?)\s+"  # Element symbol (Tag)
        r"[\d.]+\s+"  # Charge (ignored)
        r"([-\d.]+)\s+"  # X coordinate
        r"([-\d.]+)\s+"  # Y coordinate
        r"([-\d.]+)",  # Z coordinate
        re.MULTILINE,
    )

    for line_match in atom_line_pattern.finditer(coord_section):
        element = line_match.group(1)
        x = float(line_match.group(2))
        y = float(line_match.group(3))
        z = float(line_match.group(4))
        xyz_lines.append(f"{element:2s} {x:14.8f} {y:14.8f} {z:14.8f}")

    if not xyz_lines:
        raise RuntimeError(
            "Found geometry section but could not parse atom coordinates. "
            "Check NWChem output format."
        )

    return "\n".join(xyz_lines)


def parse_shielding_output(output_text: str) -> dict:
    """Parse NMR shielding tensors from NWChem output.

    Extracts isotropic chemical shielding values for each atom from the
    GIAO (Gauge-Including Atomic Orbitals) NMR calculation output.

    Args:
        output_text: Full NWChem output file content as string.

    Returns:
        Dictionary with keys matching shifts.py expected format:
        {
            "index": [1, 2, 3, ...],      # 1-based atom indices
            "atom": ["C", "H", "H", ...], # Element symbols
            "shielding": [183.4, 29.1, ...]  # Isotropic shielding in ppm
        }

    Raises:
        RuntimeError: If no shielding tensor section is found.

    Example:
        >>> output = Path("shielding.out").read_text()
        >>> data = parse_shielding_output(output)
        >>> data["atom"]  # ['C', 'H', 'H', 'H', 'H']
        >>> data["shielding"]  # [183.4, 29.1, 29.1, 29.1, 29.1]
    """
    # Check for shielding section
    # NWChem uses "GIAO Chemical Shielding Tensors" or similar header
    if "Chemical Shielding" not in output_text and "chemical shielding" not in output_text.lower():
        raise RuntimeError(
            "Could not find NMR shielding section in NWChem output. "
            "Expected 'GIAO Chemical Shielding Tensors' section. "
            "Ensure the calculation used 'property; shielding; end' and 'task dft property'."
        )

    # Initialize result lists
    indices = []
    atoms = []
    shieldings = []

    # Pattern to match atom blocks and extract isotropic shielding
    # NWChem format:
    #   Atom:    1  C
    #      isotropic  =     183.4567
    #
    # Match "Atom:" followed by index, element, then find isotropic value
    atom_block_pattern = re.compile(
        r"Atom:\s+(\d+)\s+([A-Za-z][a-z]?)"  # "Atom:    1  C"
        r".*?"  # Any content between
        r"isotropic\s*=\s*([-\d.]+)",  # "isotropic  =     183.4567"
        re.DOTALL | re.IGNORECASE,
    )

    for match in atom_block_pattern.finditer(output_text):
        atom_index = int(match.group(1))
        element = match.group(2)
        isotropic = float(match.group(3))

        indices.append(atom_index)
        atoms.append(element)
        shieldings.append(isotropic)

    if not indices:
        # Try alternative pattern for different NWChem output formats
        # Some versions might format differently
        alt_pattern = re.compile(
            r"(\d+)\s+([A-Za-z][a-z]?)\s+.*?"
            r"(?:isotropic|iso)\s*[=:]\s*([-\d.]+)",
            re.IGNORECASE,
        )

        for match in alt_pattern.finditer(output_text):
            atom_index = int(match.group(1))
            element = match.group(2)
            isotropic = float(match.group(3))

            indices.append(atom_index)
            atoms.append(element)
            shieldings.append(isotropic)

    if not indices:
        raise RuntimeError(
            "Found shielding section but could not parse atom shielding values. "
            "Expected format: 'Atom: N Element' followed by 'isotropic = value'."
        )

    return {
        "index": indices,
        "atom": atoms,
        "shielding": shieldings,
    }
