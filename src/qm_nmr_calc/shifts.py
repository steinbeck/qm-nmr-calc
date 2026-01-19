"""Shielding-to-shift conversion for NMR calculations."""

# Scaling factors for B3LYP/6-311+G(2d,p) with TMS reference.
# Reference: Pierens, J. Comput. Chem. 2014, 35, 1388-1394
# Linear regression: shift = m * shielding + b
SCALING_FACTORS: dict[str, dict[str, float]] = {
    "H": {"m": -1.0, "b": 31.8},  # TMS reference ~31.8 ppm shielding for 1H
    "C": {"m": -1.0, "b": 182.5},  # TMS reference ~182.5 ppm shielding for 13C
}


def shielding_to_shift(
    shielding_data: dict, scaling: dict[str, dict[str, float]] = SCALING_FACTORS
) -> dict[str, list[dict]]:
    """Convert shielding values to chemical shifts.

    Uses linear regression (shift = m * shielding + b) with TMS as reference.

    Args:
        shielding_data: Dict with keys 'index', 'atom', 'shielding' from NWChemParser.
            Format: {'index': [1,2,3...], 'atom': ['H','C'...], 'shielding': [29.1, 156.2...]}
        scaling: Dict mapping atom type to {'m': slope, 'b': intercept}.
            Defaults to SCALING_FACTORS for B3LYP/6-311+G(2d,p).

    Returns:
        Dict with '1H' and '13C' keys, each containing a list of shift dicts
        sorted by shift value descending (highest ppm first, standard NMR convention).
        Each shift dict: {'index': int, 'atom': str, 'shielding': float, 'shift': float}

    Raises:
        ValueError: If shielding_data is missing required keys.
    """
    # Validate input
    required_keys = {"index", "atom", "shielding"}
    missing_keys = required_keys - set(shielding_data.keys())
    if missing_keys:
        raise ValueError(
            f"shielding_data missing required keys: {', '.join(sorted(missing_keys))}"
        )

    shifts: list[dict] = []
    for idx, atom, shield in zip(
        shielding_data["index"],
        shielding_data["atom"],
        shielding_data["shielding"],
    ):
        if atom in scaling:
            m = scaling[atom]["m"]
            b = scaling[atom]["b"]
            shift = m * shield + b
            shifts.append(
                {
                    "index": idx,
                    "atom": atom,
                    "shielding": shield,
                    "shift": round(shift, 2),
                }
            )

    # Separate by nucleus type
    h_shifts = [s for s in shifts if s["atom"] == "H"]
    c_shifts = [s for s in shifts if s["atom"] == "C"]

    # Sort by shift descending (highest ppm first, standard NMR convention)
    return {
        "1H": sorted(h_shifts, key=lambda x: x["shift"], reverse=True),
        "13C": sorted(c_shifts, key=lambda x: x["shift"], reverse=True),
    }
