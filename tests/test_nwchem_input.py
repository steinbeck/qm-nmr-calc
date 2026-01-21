"""Tests for NWChem input file generation."""

import pytest

from qm_nmr_calc.nwchem import generate_optimization_input, generate_shielding_input


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

    def test_optimization_input_chcl3_dielectric(self):
        """Verify dielec 4.8 for chcl3 solvent."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
        )

        assert "dielec 4.8" in result

    def test_optimization_input_dmso_dielectric(self):
        """Verify dielec 46.0 for dmso solvent."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="dmso",
        )

        assert "dielec 46.0" in result

    def test_optimization_input_quoted_basis_set(self):
        """Verify basis set name is quoted (important for special chars like *)."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
        )

        assert 'library "6-31G*"' in result

    def test_optimization_input_max_iter(self):
        """Verify iterations parameter appears in dft block."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
            max_iter=200,
        )

        assert "iterations 200" in result

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
        assert "dielec 46.0" in result


class TestInvalidSolvent:
    """Tests for invalid solvent handling."""

    def test_invalid_solvent_raises_valueerror(self):
        """Verify ValueError with helpful message."""
        with pytest.raises(ValueError) as exc_info:
            generate_optimization_input(
                geometry_xyz="C 0.0 0.0 0.0",
                functional="b3lyp",
                basis_set="6-31G*",
                solvent="water",
            )

        error_message = str(exc_info.value).lower()
        assert "water" in error_message
        assert "chcl3" in error_message
        assert "dmso" in error_message

    def test_invalid_solvent_raises_for_shielding(self):
        """Verify ValueError for shielding input as well."""
        with pytest.raises(ValueError) as exc_info:
            generate_shielding_input(
                geometry_xyz="C 0.0 0.0 0.0",
                functional="b3lyp",
                basis_set="6-311+G(2d,p)",
                solvent="acetone",
            )

        error_message = str(exc_info.value).lower()
        assert "acetone" in error_message
        assert "chcl3" in error_message
        assert "dmso" in error_message

    def test_solvent_case_insensitive(self):
        """Both 'CHCL3' and 'chcl3' should work."""
        # Uppercase
        result_upper = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="CHCL3",
        )
        assert "dielec 4.8" in result_upper

        # Mixed case
        result_mixed = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="ChCl3",
        )
        assert "dielec 4.8" in result_mixed

        # Lowercase (standard)
        result_lower = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
        )
        assert "dielec 4.8" in result_lower
