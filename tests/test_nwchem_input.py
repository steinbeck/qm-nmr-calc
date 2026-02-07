"""Tests for NWChem input file generation."""

import pytest

from qm_nmr_calc.nwchem import generate_optimization_input, generate_shielding_input
from qm_nmr_calc.nwchem.input_gen import SUPPORTED_SOLVENTS


class TestGenerateOptimizationInput:
    """Tests for generate_optimization_input function."""

    def test_optimization_input_has_required_directives(self):
        """Check for geometry, basis, dft, cosmo, task dft optimize."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0\nH 1.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
        )

        assert "geometry" in result.lower()
        assert "basis" in result.lower()
        assert "dft" in result.lower()
        assert "cosmo" in result.lower()
        assert "task dft optimize" in result

    def test_optimization_input_chcl3_solvent(self):
        """Verify COSMO block with chcl3 solvent name."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
        )

        assert "cosmo" in result.lower()
        assert "solvent chcl3" in result.lower()

    def test_optimization_input_dmso_solvent(self):
        """Verify COSMO block with dmso solvent name."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="dmso",
        )

        assert "cosmo" in result.lower()
        assert "solvent dmso" in result.lower()

    def test_optimization_input_basis_set_in_library(self):
        """Verify basis set name appears in library directive."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
        )

        assert "library 6-31G*" in result

    def test_optimization_input_max_iter(self):
        """Verify maxiter parameter appears in driver block."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
            max_iter=200,
        )

        assert "maxiter 200" in result

    def test_optimization_input_noautosym(self):
        """Verify geometry block has noautosym."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
        )

        assert "noautosym" in result


class TestGenerateShieldingInput:
    """Tests for generate_shielding_input function."""

    def test_shielding_input_has_property_block(self):
        """Check for property/shielding block."""
        result = generate_shielding_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-311+G(2d,p)",
            solvent="chcl3",
        )

        assert "property" in result.lower()
        assert "shielding" in result.lower()

    def test_shielding_input_task_dft_property(self):
        """Check for task dft property."""
        result = generate_shielding_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-311+G(2d,p)",
            solvent="chcl3",
        )

        assert "task dft property" in result

    def test_shielding_input_has_cosmo(self):
        """Verify COSMO block present (both steps need it)."""
        result = generate_shielding_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-311+G(2d,p)",
            solvent="dmso",
        )

        assert "cosmo" in result.lower()
        assert "solvent dmso" in result.lower()


class TestVacuumSupport:
    """Tests for vacuum (gas-phase) calculations without COSMO."""

    def test_vacuum_optimization_no_cosmo(self):
        """Verify vacuum optimization has no COSMO block."""
        result = generate_optimization_input(
            geometry_xyz="H 0.0 0.0 0.0\nH 0.0 0.0 1.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="vacuum",
        )

        assert "cosmo" not in result.lower()
        assert "task dft optimize" in result

    def test_vacuum_shielding_no_cosmo(self):
        """Verify vacuum shielding has no COSMO block."""
        result = generate_shielding_input(
            geometry_xyz="H 0.0 0.0 0.0\nH 0.0 0.0 1.0",
            functional="b3lyp",
            basis_set="6-311+G(2d,p)",
            solvent="vacuum",
        )

        assert "cosmo" not in result.lower()
        assert "task dft property" in result

    def test_vacuum_case_insensitive(self):
        """Vacuum should work with any case."""
        for vacuum_str in ["vacuum", "VACUUM", "Vacuum"]:
            result = generate_optimization_input(
                geometry_xyz="H 0.0 0.0 0.0",
                functional="b3lyp",
                basis_set="6-31G*",
                solvent=vacuum_str,
            )
            assert "cosmo" not in result.lower()


class TestInvalidSolvent:
    """Tests for invalid solvent handling."""

    def test_invalid_solvent_raises_valueerror(self):
        """Verify ValueError with helpful message for unknown solvent."""
        with pytest.raises(ValueError) as exc_info:
            generate_optimization_input(
                geometry_xyz="C 0.0 0.0 0.0",
                functional="b3lyp",
                basis_set="6-31G*",
                solvent="unknown_solvent",
            )

        error_message = str(exc_info.value).lower()
        assert "unknown_solvent" in error_message
        # Should mention valid options
        assert "chcl3" in error_message or "valid" in error_message

    def test_invalid_solvent_raises_for_shielding(self):
        """Verify ValueError for shielding input as well."""
        with pytest.raises(ValueError) as exc_info:
            generate_shielding_input(
                geometry_xyz="C 0.0 0.0 0.0",
                functional="b3lyp",
                basis_set="6-311+G(2d,p)",
                solvent="not_a_solvent",
            )

        error_message = str(exc_info.value).lower()
        assert "not_a_solvent" in error_message

    def test_solvent_case_insensitive(self):
        """Both 'CHCL3' and 'chcl3' should work."""
        # Uppercase
        result_upper = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="CHCL3",
        )
        assert "cosmo" in result_upper.lower()
        assert "solvent chcl3" in result_upper.lower()

        # Mixed case
        result_mixed = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="ChCl3",
        )
        assert "cosmo" in result_mixed.lower()
        assert "solvent chcl3" in result_mixed.lower()

        # Lowercase (standard)
        result_lower = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
        )
        assert "cosmo" in result_lower.lower()
        assert "solvent chcl3" in result_lower.lower()


class TestBenzeneSolvent:
    """Tests for benzene solvent support in NWChem COSMO."""

    def test_benzene_optimization_has_cosmo(self):
        """Generate optimization input with benzene, verify COSMO block."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="benzene",
        )

        assert "cosmo" in result.lower()
        assert "solvent benzene" in result

    def test_benzene_shielding_has_cosmo(self):
        """Generate shielding input with benzene, verify COSMO block."""
        result = generate_shielding_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-311+G(2d,p)",
            solvent="benzene",
        )

        assert "cosmo" in result.lower()
        assert "solvent benzene" in result

    def test_benzene_case_insensitive(self):
        """Benzene should work with any case."""
        for benzene_str in ["benzene", "Benzene", "BENZENE"]:
            result = generate_optimization_input(
                geometry_xyz="C 0.0 0.0 0.0",
                functional="b3lyp",
                basis_set="6-31G*",
                solvent=benzene_str,
            )
            assert "cosmo" in result.lower()
            assert "solvent benzene" in result

    @pytest.mark.parametrize(
        "solvent",
        sorted(SUPPORTED_SOLVENTS - {"vacuum"}),
    )
    def test_all_supported_solvents_accepted(self, solvent):
        """All supported solvents should generate valid optimization input."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent=solvent,
        )

        assert "cosmo" in result.lower()
        assert f"solvent {solvent}" in result
