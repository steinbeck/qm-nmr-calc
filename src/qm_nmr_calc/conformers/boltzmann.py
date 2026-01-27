"""Boltzmann weight calculation for conformer ensembles."""

import math

from qm_nmr_calc.models import AtomShift, ConformerEnsemble, EnergyUnit, NMRResults

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


def average_nmr_shifts(
    per_conformer_shifts: list[list[AtomShift]],
    weights: list[float],
) -> list[AtomShift]:
    """
    Calculate population-weighted average chemical shifts across conformers.

    Args:
        per_conformer_shifts: List of shift lists, one per conformer
        weights: List of Boltzmann weights (must sum to 1.0)

    Returns:
        List of averaged AtomShift objects, sorted by shift descending

    Raises:
        ValueError: If lengths don't match or inputs are invalid

    Examples:
        >>> conf1 = [AtomShift(index=1, atom="H", shielding=30.0, shift=2.0)]
        >>> conf2 = [AtomShift(index=1, atom="H", shielding=32.0, shift=1.5)]
        >>> result = average_nmr_shifts([conf1, conf2], [0.5, 0.5])
        >>> abs(result[0].shift - 1.75) < 0.01
        True
    """
    # Validate inputs
    if len(per_conformer_shifts) != len(weights):
        raise ValueError("per_conformer_shifts and weights must have the same length")

    # Fast path: single conformer
    if len(per_conformer_shifts) == 1:
        # Return sorted copy (don't mutate input)
        return sorted(per_conformer_shifts[0], key=lambda x: x.shift, reverse=True)

    # Handle empty shifts case
    if all(len(shifts) == 0 for shifts in per_conformer_shifts):
        return []

    # Accumulate weighted shifts by atom index
    # {index: {"shielding": sum, "shift": sum, "atom": str}}
    accumulated: dict[int, dict[str, float | str]] = {}

    for i, conformer_shifts in enumerate(per_conformer_shifts):
        weight = weights[i]
        for atom_shift in conformer_shifts:
            if atom_shift.index not in accumulated:
                accumulated[atom_shift.index] = {
                    "shielding": 0.0,
                    "shift": 0.0,
                    "atom": atom_shift.atom,
                }
            accumulated[atom_shift.index]["shielding"] += weight * atom_shift.shielding
            accumulated[atom_shift.index]["shift"] += weight * atom_shift.shift

    # Build result list with rounded values
    result = []
    for index, data in accumulated.items():
        # Round: shielding to 4 decimals, shift to 2 decimals (matching shifts.py conventions)
        avg_shift = AtomShift(
            index=index,
            atom=str(data["atom"]),
            shielding=round(float(data["shielding"]), 4),
            shift=round(float(data["shift"]), 2),
        )
        result.append(avg_shift)

    # Sort by shift descending (NMR convention)
    return sorted(result, key=lambda x: x.shift, reverse=True)


def average_ensemble_nmr(
    ensemble: ConformerEnsemble,
    per_conformer_nmr: list[NMRResults],
) -> NMRResults:
    """
    Calculate Boltzmann-weighted average NMR shifts for a conformer ensemble.

    High-level orchestration function that:
    1. Calculates Boltzmann weights from conformer energies
    2. Populates ConformerData.weight for each conformer (mutates ensemble)
    3. Averages H1 and C13 shifts independently
    4. Returns single NMRResults with averaged shifts

    Args:
        ensemble: ConformerEnsemble with energies and temperature
        per_conformer_nmr: List of NMRResults, one per conformer

    Returns:
        NMRResults with averaged shifts and metadata from first conformer

    Raises:
        ValueError: If inputs are invalid or missing required data

    Examples:
        >>> ensemble = ConformerEnsemble(
        ...     method="rdkit_kdg",
        ...     conformers=[
        ...         ConformerData(conformer_id="c1", energy=0.0, energy_unit="kcal_mol"),
        ...         ConformerData(conformer_id="c2", energy=0.59, energy_unit="kcal_mol"),
        ...     ],
        ...     temperature_k=298.15,
        ... )
        >>> nmr1 = NMRResults(h1_shifts=[...], c13_shifts=[...], ...)
        >>> nmr2 = NMRResults(h1_shifts=[...], c13_shifts=[...], ...)
        >>> result = average_ensemble_nmr(ensemble, [nmr1, nmr2])
        >>> abs(ensemble.conformers[0].weight - 0.73) < 0.01
        True
    """
    # Validate inputs
    if not per_conformer_nmr:
        raise ValueError("Must provide at least one conformer NMR result")
    if len(ensemble.conformers) != len(per_conformer_nmr):
        raise ValueError("Must have same number of conformers and NMR results")

    # Validate all conformers have energies
    for i, conformer in enumerate(ensemble.conformers):
        if conformer.energy is None:
            raise ValueError("All conformers must have energies")
        if conformer.energy_unit is None:
            raise ValueError("All conformers must have energy_unit")

    # Extract energies and energy_unit
    energies = [c.energy for c in ensemble.conformers]
    energy_unit = ensemble.conformers[0].energy_unit

    # Calculate Boltzmann weights
    weights = calculate_boltzmann_weights(energies, ensemble.temperature_k, energy_unit)

    # Populate weights in ensemble (mutate in place)
    for i, conformer in enumerate(ensemble.conformers):
        conformer.weight = weights[i]

    # Collect per-conformer shifts by nucleus type
    per_conformer_h1 = [nmr.h1_shifts for nmr in per_conformer_nmr]
    per_conformer_c13 = [nmr.c13_shifts for nmr in per_conformer_nmr]

    # Average shifts
    avg_h1 = average_nmr_shifts(per_conformer_h1, weights)
    avg_c13 = average_nmr_shifts(per_conformer_c13, weights)

    # Return NMRResults with averaged shifts and metadata from first conformer
    return NMRResults(
        h1_shifts=avg_h1,
        c13_shifts=avg_c13,
        functional=per_conformer_nmr[0].functional,
        basis_set=per_conformer_nmr[0].basis_set,
        solvent=per_conformer_nmr[0].solvent,
    )
