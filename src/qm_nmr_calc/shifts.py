"""Shielding-to-shift conversion for NMR calculations.

Scaling factors derived from DELTA50 benchmark dataset using linear regression
between NWChem-calculated shieldings and experimental shifts.

Formula: shift = slope * shielding + intercept
"""

from functools import cache
from importlib.resources import files

import orjson


@cache
def load_scaling_factors() -> dict:
    """Load scaling factors from package data (cached after first call).

    Returns:
        Dict mapping factor keys to factor dicts.
        Key format: "B3LYP/6-311+G(2d,p)/1H/CHCl3"
        Value: {slope, intercept, r_squared, mae, rmsd, n_points, ...}
    """
    data_dir = files("qm_nmr_calc").joinpath("data")
    json_bytes = data_dir.joinpath("scaling_factors.json").read_bytes()
    return orjson.loads(json_bytes)


def get_scaling_factor(
    functional: str,
    basis_set: str,
    nucleus: str,
    solvent: str
) -> dict:
    """Get scaling factor for specific calculation parameters.

    Args:
        functional: DFT functional (e.g., 'B3LYP')
        basis_set: Basis set (e.g., '6-311+G(2d,p)')
        nucleus: Nucleus type ('1H' or '13C')
        solvent: Solvent name (e.g., 'CHCl3', 'DMSO', 'vacuum')

    Returns:
        Dict with: slope, intercept, r_squared, mae, rmsd, n_points

    Raises:
        ValueError: If no factor exists for this combination
    """
    # Normalize solvent name to match scaling factors format
    solvent_map = {"chcl3": "CHCl3", "dmso": "DMSO", "vacuum": "vacuum"}
    normalized_solvent = solvent_map.get(solvent.lower(), solvent)

    key = f"{functional}/{basis_set}/{nucleus}/{normalized_solvent}"
    factors = load_scaling_factors()

    if key not in factors:
        available = list(factors.keys())
        raise ValueError(
            f"No scaling factor for {key}.\n"
            f"Supported combinations: {available}"
        )

    return factors[key]


def shielding_to_shift(
    shielding_data: dict,
    functional: str,
    basis_set: str,
    solvent: str,
) -> dict[str, list[dict]]:
    """Convert shielding values to chemical shifts.

    Uses empirical linear regression: shift = slope * shielding + intercept
    Scaling factors are from DELTA50 benchmark, fitted to experimental data.

    Args:
        shielding_data: Dict with keys 'index', 'atom', 'shielding' from NWChemParser.
            Format: {'index': [1,2,3...], 'atom': ['H','C'...], 'shielding': [29.1, 156.2...]}
        functional: DFT functional (e.g., 'B3LYP')
        basis_set: Basis set (e.g., '6-311+G(2d,p)')
        solvent: Solvent name (e.g., 'CHCl3', 'DMSO', 'vacuum')

    Returns:
        Dict with '1H' and '13C' keys, each containing a list of shift dicts
        sorted by shift value descending (highest ppm first, standard NMR convention).
        Each shift dict: {'index': int, 'atom': str, 'shielding': float, 'shift': float}

    Raises:
        ValueError: If shielding_data is missing required keys or no scaling factor exists.
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
        # Determine nucleus type
        nucleus = "1H" if atom == "H" else "13C" if atom == "C" else None
        if nucleus is None:
            continue  # Skip unsupported atom types

        # Get scaling factor for this nucleus/functional/basis_set/solvent
        factor = get_scaling_factor(functional, basis_set, nucleus, solvent)

        # Apply regression formula
        slope = factor["slope"]
        intercept = factor["intercept"]
        shift = slope * shield + intercept

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
