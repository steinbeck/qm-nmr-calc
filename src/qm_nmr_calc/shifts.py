"""Shielding-to-shift conversion for NMR calculations."""

# Empirical scaling factors from CHESHIRE Chemical Shift Repository
# Source: http://cheshirenmr.info/ScalingFactors.htm (Pierens 2014)
# Table #1b: Chloroform DFT methods (G09 - SMD solvation model)
#
# Formula: shift = slope * shielding + intercept
#
# Note: These are empirical linear regression parameters fitted to experimental
# data, NOT simple TMS referencing. The slope accounts for systematic errors
# in the DFT method.

SCALING_FACTORS: dict[str, dict[str, dict[str, float]]] = {
    "draft": {
        # B3LYP/6-31G(d)//B3LYP/6-31G(d) gas-phase
        # Using same basis for geometry and NMR
        "H": {"m": -1.0157, "b": 32.2109},
        "C": {"m": -0.9449, "b": 188.4418},
    },
    "production": {
        # B3LYP/6-31+G(d,p)//B3LYP/6-311+G(2d,p) gas-phase (CHESHIRE Table #1a)
        # Closest match to our setup: 6-31G* geometry, 6-311+G(2d,p) NMR
        "H": {"m": -1.0592, "b": 31.9654},
        "C": {"m": -1.0311, "b": 180.7713},
    },
}


def get_scaling_factors(preset: str = "production") -> dict[str, dict[str, float]]:
    """Get empirical scaling factors for a given preset.

    Args:
        preset: Preset name ('draft' or 'production')

    Returns:
        Dict mapping atom type to {'m': slope, 'b': intercept}
    """
    return SCALING_FACTORS.get(preset, SCALING_FACTORS["production"])


def shielding_to_shift(
    shielding_data: dict,
    scaling: dict[str, dict[str, float]] = None,
    preset: str = "production",
) -> dict[str, list[dict]]:
    """Convert shielding values to chemical shifts.

    Uses empirical linear regression: shift = m * shielding + b
    Scaling factors are from CHESHIRE (Pierens 2014), fitted to experimental data.

    Args:
        shielding_data: Dict with keys 'index', 'atom', 'shielding' from NWChemParser.
            Format: {'index': [1,2,3...], 'atom': ['H','C'...], 'shielding': [29.1, 156.2...]}
        scaling: Optional custom scaling factors. If None, uses preset-specific values.
        preset: Preset name ('draft' or 'production') to select scaling factors.

    Returns:
        Dict with '1H' and '13C' keys, each containing a list of shift dicts
        sorted by shift value descending (highest ppm first, standard NMR convention).
        Each shift dict: {'index': int, 'atom': str, 'shielding': float, 'shift': float}

    Raises:
        ValueError: If shielding_data is missing required keys.
    """
    # Use preset-specific empirical scaling factors if no custom scaling provided
    if scaling is None:
        scaling = get_scaling_factors(preset)

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
