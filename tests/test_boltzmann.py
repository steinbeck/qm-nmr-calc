"""Tests for Boltzmann weight calculation."""

import math

import pytest

from qm_nmr_calc.conformers.boltzmann import (
    average_ensemble_nmr,
    average_nmr_shifts,
    calculate_boltzmann_weights,
)
from qm_nmr_calc.models import AtomShift, ConformerData, ConformerEnsemble, NMRResults


class TestBasicBehavior:
    """Test basic Boltzmann weight calculation behavior."""

    def test_single_conformer_returns_weight_one(self):
        """Single conformer should have weight = 1.0."""
        weights = calculate_boltzmann_weights([0.0])
        assert weights == [1.0]

        weights = calculate_boltzmann_weights([42.5])
        assert weights == [1.0]

        weights = calculate_boltzmann_weights([-100.0])
        assert weights == [1.0]

    def test_equal_energies_return_equal_weights(self):
        """Equal energies should produce equal weights (1/N each)."""
        weights = calculate_boltzmann_weights([5.0, 5.0, 5.0])
        expected = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
        assert len(weights) == 3
        for w, exp in zip(weights, expected):
            assert abs(w - exp) < 1e-10

    def test_two_conformers_known_case(self):
        """Two conformers 0.59 kcal/mol apart at 298.15 K -> 73%/27%."""
        # RT = 0.001987204 * 298.15 = 0.5922 kcal/mol
        # exp(-0.59/0.5922) = exp(-0.9966) = 0.3691
        # weights = [1.0, 0.3691] / (1.0 + 0.3691) = [0.7305, 0.2695]
        weights = calculate_boltzmann_weights([0.0, 0.59])
        assert len(weights) == 2
        assert abs(weights[0] - 0.7305) < 0.001
        assert abs(weights[1] - 0.2695) < 0.001
        assert abs(sum(weights) - 1.0) < 1e-10

    def test_large_energy_difference_dominant_conformer(self):
        """Large energy difference (20 kcal/mol) -> dominant conformer weight near 1.0."""
        weights = calculate_boltzmann_weights([0.0, 20.0])
        assert len(weights) == 2
        assert weights[0] > 0.99999  # Essentially 1.0
        assert weights[1] < 0.00001  # Essentially 0.0
        assert abs(sum(weights) - 1.0) < 1e-10

    def test_five_conformers_monotonic_decreasing(self):
        """Five conformers with increasing energies should have monotonically decreasing weights."""
        energies = [0.0, 0.5, 1.0, 2.0, 5.0]
        weights = calculate_boltzmann_weights(energies)
        assert len(weights) == 5
        # Check monotonic decrease
        for i in range(len(weights) - 1):
            assert weights[i] > weights[i + 1]
        # Check sum
        assert abs(sum(weights) - 1.0) < 1e-10


class TestEnergyUnitConversion:
    """Test energy unit conversion handling."""

    def test_kcal_mol_no_conversion(self):
        """kcal_mol units should be used directly."""
        weights1 = calculate_boltzmann_weights([10.0, 11.0], energy_unit="kcal_mol")
        weights2 = calculate_boltzmann_weights([0.0, 1.0], energy_unit="kcal_mol")
        # Same relative difference, should give same weights
        assert len(weights1) == len(weights2) == 2
        for w1, w2 in zip(weights1, weights2):
            assert abs(w1 - w2) < 1e-10

    def test_hartree_conversion(self):
        """Hartree energies should be converted to kcal/mol."""
        # 1 hartree = 627.5095 kcal/mol
        # Energy difference of 0.001 hartree = 0.6275 kcal/mol
        weights_hartree = calculate_boltzmann_weights([-76.0, -75.999], energy_unit="hartree")
        weights_kcal = calculate_boltzmann_weights([0.0, 0.6275], energy_unit="kcal_mol")

        assert len(weights_hartree) == 2
        # Should give approximately the same weights as equivalent kcal/mol difference
        assert abs(weights_hartree[0] - weights_kcal[0]) < 0.001
        assert abs(weights_hartree[1] - weights_kcal[1]) < 0.001

    def test_kj_mol_conversion(self):
        """kJ/mol energies should be converted to kcal/mol."""
        # 1 kcal = 4.184 kJ
        # 0.59 kcal/mol = 2.468 kJ/mol
        weights_kj = calculate_boltzmann_weights([41.84, 44.31], energy_unit="kj_mol")
        weights_kcal = calculate_boltzmann_weights([10.0, 10.59], energy_unit="kcal_mol")

        assert len(weights_kj) == 2
        # Should give approximately the same weights as equivalent kcal/mol
        assert abs(weights_kj[0] - weights_kcal[0]) < 0.001
        assert abs(weights_kj[1] - weights_kcal[1]) < 0.001


class TestNumericalStability:
    """Test numerical stability with extreme values."""

    def test_extreme_energy_range_no_overflow(self):
        """Extreme energy range [0.0, 100.0] should not overflow."""
        weights = calculate_boltzmann_weights([0.0, 100.0])
        assert len(weights) == 2
        assert all(math.isfinite(w) for w in weights)
        assert abs(sum(weights) - 1.0) < 1e-10
        assert weights[0] > 0.999  # First dominates

    def test_large_negative_energies_relative_only(self):
        """Large negative energies should give same result as relative values."""
        weights1 = calculate_boltzmann_weights([-1000.0, -999.0])
        weights2 = calculate_boltzmann_weights([0.0, 1.0])

        assert len(weights1) == len(weights2) == 2
        # Relative differences matter, not absolute values
        for w1, w2 in zip(weights1, weights2):
            assert abs(w1 - w2) < 1e-10

    def test_large_absolute_values_correct_relative_weights(self):
        """Large absolute values should still produce correct relative weights."""
        # Same relative differences as [0.0, 0.5, 1.0]
        weights1 = calculate_boltzmann_weights([500.0, 500.5, 501.0])
        weights2 = calculate_boltzmann_weights([0.0, 0.5, 1.0])

        assert len(weights1) == len(weights2) == 3
        for w1, w2 in zip(weights1, weights2):
            assert abs(w1 - w2) < 1e-10


class TestTemperature:
    """Test temperature parameter handling."""

    def test_default_temperature_298_15K(self):
        """Default temperature should be 298.15 K."""
        # This is implicitly tested by the known case at 298.15 K
        weights = calculate_boltzmann_weights([0.0, 0.59])
        # If temperature is correct, we get the expected weights
        assert abs(weights[0] - 0.7305) < 0.001

    def test_higher_temperature_more_equal_distribution(self):
        """Higher temperature should give more equal weight distribution."""
        weights_low_temp = calculate_boltzmann_weights([0.0, 1.0], temperature_k=100.0)
        weights_high_temp = calculate_boltzmann_weights([0.0, 1.0], temperature_k=500.0)

        # At higher temperature, the second conformer gets more weight
        assert weights_high_temp[1] > weights_low_temp[1]
        # But first still dominates (lower energy)
        assert weights_high_temp[0] > weights_high_temp[1]


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_list_raises_error(self):
        """Empty energy list should raise ValueError."""
        with pytest.raises(ValueError, match="At least one energy value required"):
            calculate_boltzmann_weights([])

    def test_zero_temperature_raises_error(self):
        """Zero temperature should raise ValueError."""
        with pytest.raises(ValueError, match="Temperature must be positive"):
            calculate_boltzmann_weights([0.0, 1.0], temperature_k=0.0)

    def test_negative_temperature_raises_error(self):
        """Negative temperature should raise ValueError."""
        with pytest.raises(ValueError, match="Temperature must be positive"):
            calculate_boltzmann_weights([0.0, 1.0], temperature_k=-10.0)


class TestReturnContract:
    """Test return value contract guarantees."""

    def test_weights_sum_to_one(self):
        """Weights should always sum to 1.0."""
        test_cases = [
            [0.0, 1.0, 2.0],
            [5.5, 6.2, 7.1, 8.3],
            [-10.0, -9.0, -8.5],
            [100.0, 101.0, 102.0, 103.0, 104.0],
        ]
        for energies in test_cases:
            weights = calculate_boltzmann_weights(energies)
            assert abs(sum(weights) - 1.0) < 1e-10

    def test_all_weights_non_negative(self):
        """All weights should be >= 0.0."""
        test_cases = [
            [0.0, 1.0],
            [-50.0, -49.0, -48.0],
            [10.0, 15.0, 20.0, 25.0],
        ]
        for energies in test_cases:
            weights = calculate_boltzmann_weights(energies)
            assert all(w >= 0.0 for w in weights)

    def test_lower_energy_higher_weight(self):
        """Lower energy conformers should have higher weight."""
        energies = [1.0, 2.0, 3.0, 4.0, 5.0]
        weights = calculate_boltzmann_weights(energies)
        # Sorted energies should give sorted weights (descending)
        for i in range(len(weights) - 1):
            assert weights[i] > weights[i + 1]

    def test_return_length_matches_input(self):
        """Returned list should have same length as input."""
        for n in [1, 2, 5, 10, 20]:
            energies = list(range(n))
            weights = calculate_boltzmann_weights(energies)
            assert len(weights) == n


class TestAverageNMRShifts:
    """Test average_nmr_shifts function for shift averaging."""

    def test_equal_weights_arithmetic_mean(self):
        """Two conformers with equal weights [0.5, 0.5] should return arithmetic mean."""
        conf1 = [AtomShift(index=1, atom="H", shielding=30.0, shift=2.0)]
        conf2 = [AtomShift(index=1, atom="H", shielding=32.0, shift=1.5)]
        weights = [0.5, 0.5]

        result = average_nmr_shifts([conf1, conf2], weights)

        assert len(result) == 1
        assert result[0].index == 1
        assert result[0].atom == "H"
        assert abs(result[0].shielding - 31.0) < 1e-4
        assert abs(result[0].shift - 1.75) < 0.01

    def test_unequal_weights_weighted_average(self):
        """Two conformers with unequal weights [0.75, 0.25] should return weighted average."""
        # shielding: 0.75*30 + 0.25*34 = 22.5 + 8.5 = 31.0
        # shift: 0.75*2.0 + 0.25*1.0 = 1.5 + 0.25 = 1.75
        conf1 = [AtomShift(index=1, atom="H", shielding=30.0, shift=2.0)]
        conf2 = [AtomShift(index=1, atom="H", shielding=34.0, shift=1.0)]
        weights = [0.75, 0.25]

        result = average_nmr_shifts([conf1, conf2], weights)

        assert len(result) == 1
        assert abs(result[0].shielding - 31.0) < 1e-4
        assert abs(result[0].shift - 1.75) < 0.01

    def test_single_conformer_identity(self):
        """Single conformer with weight 1.0 should return identical shifts."""
        conf = [
            AtomShift(index=1, atom="H", shielding=30.5, shift=2.25),
            AtomShift(index=2, atom="C", shielding=150.0, shift=45.0),
        ]
        weights = [1.0]

        result = average_nmr_shifts([conf], weights)

        assert len(result) == 2
        # Results sorted by shift descending, so C (45.0) before H (2.25)
        assert result[0].index == 2
        assert result[0].atom == "C"
        assert abs(result[0].shielding - 150.0) < 1e-4
        assert abs(result[0].shift - 45.0) < 0.01
        assert result[1].index == 1
        assert result[1].atom == "H"
        assert abs(result[1].shielding - 30.5) < 1e-4
        assert abs(result[1].shift - 2.25) < 0.01

    def test_multiple_atoms_averaged_independently(self):
        """Multiple atoms (H and C) should be averaged independently."""
        conf1 = [
            AtomShift(index=1, atom="H", shielding=30.0, shift=2.0),
            AtomShift(index=2, atom="H", shielding=28.0, shift=4.0),
            AtomShift(index=3, atom="C", shielding=150.0, shift=50.0),
        ]
        conf2 = [
            AtomShift(index=1, atom="H", shielding=32.0, shift=1.0),
            AtomShift(index=2, atom="H", shielding=30.0, shift=3.0),
            AtomShift(index=3, atom="C", shielding=160.0, shift=45.0),
        ]
        weights = [0.6, 0.4]

        result = average_nmr_shifts([conf1, conf2], weights)

        assert len(result) == 3
        # Results sorted by shift descending: C(48.0), H2(3.6), H1(1.6)
        # C atom 3: shielding = 0.6*150 + 0.4*160 = 154.0, shift = 0.6*50 + 0.4*45 = 48.0
        assert result[0].index == 3
        assert result[0].atom == "C"
        assert abs(result[0].shielding - 154.0) < 1e-4
        assert abs(result[0].shift - 48.0) < 0.01
        # H atom 2: shielding = 0.6*28 + 0.4*30 = 28.8, shift = 0.6*4.0 + 0.4*3.0 = 3.6
        assert result[1].index == 2
        assert result[1].atom == "H"
        assert abs(result[1].shielding - 28.8) < 1e-4
        assert abs(result[1].shift - 3.6) < 0.01
        # H atom 1: shielding = 0.6*30 + 0.4*32 = 30.8, shift = 0.6*2.0 + 0.4*1.0 = 1.6
        assert result[2].index == 1
        assert result[2].atom == "H"
        assert abs(result[2].shielding - 30.8) < 1e-4
        assert abs(result[2].shift - 1.6) < 0.01

    def test_mismatched_count_raises_error(self):
        """Mismatched conformer count and weight count should raise ValueError."""
        conf1 = [AtomShift(index=1, atom="H", shielding=30.0, shift=2.0)]
        conf2 = [AtomShift(index=1, atom="H", shielding=32.0, shift=1.5)]
        weights = [1.0]  # Only one weight for two conformers

        with pytest.raises(ValueError, match="per_conformer_shifts.*weights"):
            average_nmr_shifts([conf1, conf2], weights)

    def test_empty_shifts_returns_empty(self):
        """Empty shift lists should return empty result."""
        result = average_nmr_shifts([[], []], [0.5, 0.5])
        assert result == []

    def test_sorted_descending_by_shift(self):
        """Result should be sorted by shift descending (NMR convention)."""
        conf1 = [
            AtomShift(index=1, atom="H", shielding=30.0, shift=1.0),
            AtomShift(index=2, atom="H", shielding=28.0, shift=3.0),
            AtomShift(index=3, atom="H", shielding=29.0, shift=2.0),
        ]
        conf2 = [
            AtomShift(index=1, atom="H", shielding=31.0, shift=1.5),
            AtomShift(index=2, atom="H", shielding=27.0, shift=3.5),
            AtomShift(index=3, atom="H", shielding=28.0, shift=2.5),
        ]
        weights = [0.5, 0.5]

        result = average_nmr_shifts([conf1, conf2], weights)

        assert len(result) == 3
        # Should be sorted by shift descending
        assert result[0].shift > result[1].shift > result[2].shift
        # Verify it's index 2, 3, 1
        assert result[0].index == 2
        assert result[1].index == 3
        assert result[2].index == 1


class TestAverageEnsembleNMR:
    """Test average_ensemble_nmr function for high-level orchestration."""

    def test_two_conformers_known_energies(self):
        """Two conformers with known energies should produce correct weights and averages."""
        # Energies [0.0, 0.59] kcal/mol at 298.15 K -> weights [0.73, 0.27]
        ensemble = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[
                ConformerData(conformer_id="conf_001", energy=0.0, energy_unit="kcal_mol"),
                ConformerData(conformer_id="conf_002", energy=0.59, energy_unit="kcal_mol"),
            ],
            temperature_k=298.15,
        )
        nmr1 = NMRResults(
            h1_shifts=[AtomShift(index=1, atom="H", shielding=30.0, shift=2.0)],
            c13_shifts=[AtomShift(index=2, atom="C", shielding=150.0, shift=50.0)],
            functional="B3LYP",
            basis_set="6-31G*",
            solvent="chloroform",
        )
        nmr2 = NMRResults(
            h1_shifts=[AtomShift(index=1, atom="H", shielding=32.0, shift=1.0)],
            c13_shifts=[AtomShift(index=2, atom="C", shielding=160.0, shift=40.0)],
            functional="B3LYP",
            basis_set="6-31G*",
            solvent="chloroform",
        )

        result = average_ensemble_nmr(ensemble, [nmr1, nmr2])

        # Verify weights populated
        assert ensemble.conformers[0].weight is not None
        assert ensemble.conformers[1].weight is not None
        assert abs(ensemble.conformers[0].weight - 0.7305) < 0.001
        assert abs(ensemble.conformers[1].weight - 0.2695) < 0.001

        # Verify averaged shifts
        # H1: 0.73*2.0 + 0.27*1.0 = 1.46 + 0.27 = 1.73
        assert len(result.h1_shifts) == 1
        assert abs(result.h1_shifts[0].shift - 1.73) < 0.05
        # C13: 0.73*50.0 + 0.27*40.0 = 36.5 + 10.8 = 47.3
        assert len(result.c13_shifts) == 1
        assert abs(result.c13_shifts[0].shift - 47.3) < 0.5

        # Verify metadata copied
        assert result.functional == "B3LYP"
        assert result.basis_set == "6-31G*"
        assert result.solvent == "chloroform"

    def test_single_conformer_identity(self):
        """Single conformer ensemble should return identical NMRResults with weight=1.0."""
        ensemble = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[
                ConformerData(conformer_id="conf_001", energy=0.0, energy_unit="kcal_mol"),
            ],
            temperature_k=298.15,
        )
        nmr = NMRResults(
            h1_shifts=[AtomShift(index=1, atom="H", shielding=30.5, shift=2.25)],
            c13_shifts=[AtomShift(index=2, atom="C", shielding=150.0, shift=45.0)],
            functional="wB97X-D",
            basis_set="6-311+G(2d,p)",
            solvent="dmso",
        )

        result = average_ensemble_nmr(ensemble, [nmr])

        # Weight should be 1.0
        assert ensemble.conformers[0].weight == 1.0

        # Shifts should be identical
        assert len(result.h1_shifts) == 1
        assert abs(result.h1_shifts[0].shift - 2.25) < 0.01
        assert len(result.c13_shifts) == 1
        assert abs(result.c13_shifts[0].shift - 45.0) < 0.01

        # Metadata preserved
        assert result.functional == "wB97X-D"
        assert result.basis_set == "6-311+G(2d,p)"
        assert result.solvent == "dmso"

    def test_shifts_sorted_descending(self):
        """Returned shifts should be sorted descending by shift value."""
        ensemble = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[
                ConformerData(conformer_id="conf_001", energy=0.0, energy_unit="kcal_mol"),
            ],
            temperature_k=298.15,
        )
        nmr = NMRResults(
            h1_shifts=[
                AtomShift(index=1, atom="H", shielding=30.0, shift=1.0),
                AtomShift(index=2, atom="H", shielding=28.0, shift=3.0),
                AtomShift(index=3, atom="H", shielding=29.0, shift=2.0),
            ],
            c13_shifts=[
                AtomShift(index=4, atom="C", shielding=150.0, shift=40.0),
                AtomShift(index=5, atom="C", shielding=140.0, shift=60.0),
            ],
            functional="B3LYP",
            basis_set="6-31G*",
            solvent="chloroform",
        )

        result = average_ensemble_nmr(ensemble, [nmr])

        # H1 shifts sorted descending
        assert result.h1_shifts[0].shift > result.h1_shifts[1].shift > result.h1_shifts[2].shift
        assert result.h1_shifts[0].index == 2  # shift=3.0
        assert result.h1_shifts[1].index == 3  # shift=2.0
        assert result.h1_shifts[2].index == 1  # shift=1.0

        # C13 shifts sorted descending
        assert result.c13_shifts[0].shift > result.c13_shifts[1].shift
        assert result.c13_shifts[0].index == 5  # shift=60.0
        assert result.c13_shifts[1].index == 4  # shift=40.0

    def test_empty_nmr_list_raises_error(self):
        """Empty per_conformer_nmr list should raise ValueError."""
        ensemble = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[],
            temperature_k=298.15,
        )

        with pytest.raises(ValueError, match="at least one conformer"):
            average_ensemble_nmr(ensemble, [])

    def test_mismatched_count_raises_error(self):
        """Mismatched ensemble.conformers and per_conformer_nmr lengths should raise ValueError."""
        ensemble = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[
                ConformerData(conformer_id="conf_001", energy=0.0, energy_unit="kcal_mol"),
                ConformerData(conformer_id="conf_002", energy=0.5, energy_unit="kcal_mol"),
            ],
            temperature_k=298.15,
        )
        nmr = NMRResults(
            h1_shifts=[],
            c13_shifts=[],
            functional="B3LYP",
            basis_set="6-31G*",
            solvent="chloroform",
        )

        with pytest.raises(ValueError, match="same number of conformers"):
            average_ensemble_nmr(ensemble, [nmr])  # Only one NMR for two conformers

    def test_missing_energy_raises_error(self):
        """Conformers with None energy should raise ValueError."""
        ensemble = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[
                ConformerData(conformer_id="conf_001", energy=0.0, energy_unit="kcal_mol"),
                ConformerData(conformer_id="conf_002", energy=None, energy_unit="kcal_mol"),
            ],
            temperature_k=298.15,
        )
        nmr = NMRResults(
            h1_shifts=[],
            c13_shifts=[],
            functional="B3LYP",
            basis_set="6-31G*",
            solvent="chloroform",
        )

        with pytest.raises(ValueError, match="All conformers must have energies"):
            average_ensemble_nmr(ensemble, [nmr, nmr])

    def test_missing_energy_unit_raises_error(self):
        """Conformers with None energy_unit should raise ValueError."""
        ensemble = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[
                ConformerData(conformer_id="conf_001", energy=0.0, energy_unit="kcal_mol"),
                ConformerData(conformer_id="conf_002", energy=0.5, energy_unit=None),
            ],
            temperature_k=298.15,
        )
        nmr = NMRResults(
            h1_shifts=[],
            c13_shifts=[],
            functional="B3LYP",
            basis_set="6-31G*",
            solvent="chloroform",
        )

        with pytest.raises(ValueError, match="All conformers must have energy_unit"):
            average_ensemble_nmr(ensemble, [nmr, nmr])
