"""Tests for Boltzmann weight calculation."""

import math

import pytest

from qm_nmr_calc.conformers.boltzmann import calculate_boltzmann_weights


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
