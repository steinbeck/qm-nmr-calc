"""Boltzmann weight calculation for conformer ensembles."""

import math

from qm_nmr_calc.models import EnergyUnit

# Physical constants
R_KCAL = 0.001987204  # Gas constant in kcal/(mol*K)
HARTREE_TO_KCAL = 627.5095  # Conversion factor: 1 hartree = 627.5095 kcal/mol
KJ_TO_KCAL = 1.0 / 4.184  # Conversion factor: 1 kcal = 4.184 kJ


def calculate_boltzmann_weights(
    energies: list[float],
    temperature_k: float = 298.15,
    energy_unit: EnergyUnit = "kcal_mol",
) -> list[float]:
    """
    Calculate Boltzmann weights from conformer energies.

    Uses the exp-normalize trick for numerical stability:
    - Subtract minimum energy before exponentiation (all relative energies >= 0)
    - This prevents overflow (exp of large positive) and handles underflow gracefully

    Args:
        energies: List of conformer energies
        temperature_k: Temperature in Kelvin (default: 298.15 K)
        energy_unit: Unit of input energies ("kcal_mol", "hartree", or "kj_mol")

    Returns:
        List of Boltzmann weights (same length as energies), summing to 1.0

    Raises:
        ValueError: If energies list is empty or temperature is non-positive

    Examples:
        >>> calculate_boltzmann_weights([0.0])
        [1.0]
        >>> weights = calculate_boltzmann_weights([0.0, 0.59])
        >>> abs(weights[0] - 0.7305) < 0.001
        True
    """
    # Validate inputs
    if not energies:
        raise ValueError("At least one energy value required")
    if temperature_k <= 0:
        raise ValueError("Temperature must be positive")

    # Single conformer shortcut
    if len(energies) == 1:
        return [1.0]

    # Convert energies to kcal/mol
    if energy_unit == "kcal_mol":
        energies_kcal = energies
    elif energy_unit == "hartree":
        energies_kcal = [e * HARTREE_TO_KCAL for e in energies]
    elif energy_unit == "kj_mol":
        energies_kcal = [e * KJ_TO_KCAL for e in energies]
    else:
        raise ValueError(f"Unknown energy unit: {energy_unit}")

    # Calculate RT in kcal/mol
    rt = R_KCAL * temperature_k

    # Apply exp-normalize trick for numerical stability
    # 1. Find minimum energy (will have relative energy = 0, weight = exp(0) = 1.0)
    min_energy = min(energies_kcal)

    # 2. Calculate relative energies (all >= 0)
    relative_energies = [e - min_energy for e in energies_kcal]

    # 3. Calculate unnormalized weights: w_i = exp(-E_rel_i / RT)
    # All exponent arguments are <= 0, so no overflow risk
    unnormalized_weights = [math.exp(-e_rel / rt) for e_rel in relative_energies]

    # 4. Normalize so weights sum to 1.0
    total_weight = sum(unnormalized_weights)
    weights = [w / total_weight for w in unnormalized_weights]

    return weights
