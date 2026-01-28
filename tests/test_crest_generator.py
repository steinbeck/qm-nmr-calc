"""Tests for CREST conformer generation utilities."""

import unittest.mock as mock
from pathlib import Path

import pytest

from qm_nmr_calc.conformers.crest_generator import (
    ALPB_SOLVENT_MAP,
    CRESTConformer,
    detect_crest_available,
    get_alpb_solvent,
    parse_crest_ensemble,
)


# Tests for detect_crest_available()
class TestDetectCRESTAvailable:
    """Test CREST/xTB binary detection."""

    def test_both_binaries_present(self):
        """detect_crest_available returns True when both crest and xtb are found."""
        with mock.patch("shutil.which") as mock_which:
            # Both binaries found
            mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ["crest", "xtb"] else None
            result = detect_crest_available()
            assert result is True

    def test_only_crest_present(self):
        """detect_crest_available returns False when only crest is found."""
        with mock.patch("shutil.which") as mock_which:
            # Only crest found
            mock_which.side_effect = lambda cmd: "/usr/bin/crest" if cmd == "crest" else None
            result = detect_crest_available()
            assert result is False

    def test_only_xtb_present(self):
        """detect_crest_available returns False when only xtb is found."""
        with mock.patch("shutil.which") as mock_which:
            # Only xtb found
            mock_which.side_effect = lambda cmd: "/usr/bin/xtb" if cmd == "xtb" else None
            result = detect_crest_available()
            assert result is False

    def test_neither_binary_present(self):
        """detect_crest_available returns False when neither binary is found."""
        with mock.patch("shutil.which") as mock_which:
            # Neither found
            mock_which.return_value = None
            result = detect_crest_available()
            assert result is False


# Tests for get_alpb_solvent()
class TestGetALPBSolvent:
    """Test ALPB solvent mapping."""

    def test_chcl3_uppercase(self):
        """get_alpb_solvent maps 'CHCl3' to 'chcl3'."""
        result = get_alpb_solvent("CHCl3")
        assert result == "chcl3"

    def test_chcl3_lowercase(self):
        """get_alpb_solvent maps 'chcl3' to 'chcl3'."""
        result = get_alpb_solvent("chcl3")
        assert result == "chcl3"

    def test_dmso_uppercase(self):
        """get_alpb_solvent maps 'DMSO' to 'dmso'."""
        result = get_alpb_solvent("DMSO")
        assert result == "dmso"

    def test_dmso_lowercase(self):
        """get_alpb_solvent maps 'dmso' to 'dmso'."""
        result = get_alpb_solvent("dmso")
        assert result == "dmso"

    def test_vacuum_returns_none(self):
        """get_alpb_solvent returns None for 'vacuum'."""
        result = get_alpb_solvent("vacuum")
        assert result is None

    def test_unsupported_solvent_returns_none(self):
        """get_alpb_solvent returns None for unsupported solvents like 'water'."""
        result = get_alpb_solvent("water")
        assert result is None


# Tests for parse_crest_ensemble()
class TestParseCRESTEnsemble:
    """Test CREST multi-structure XYZ parsing."""

    def test_single_conformer(self, tmp_path):
        """parse_crest_ensemble correctly parses a single conformer file."""
        # Create test XYZ file with one conformer
        xyz_file = tmp_path / "single.xyz"
        xyz_content = """3
-12.34567890
C  0.0  0.0  0.0
H  1.0  0.0  0.0
H  0.0  1.0  0.0
"""
        xyz_file.write_text(xyz_content)

        # Parse
        conformers = parse_crest_ensemble(xyz_file)

        # Verify
        assert len(conformers) == 1
        assert conformers[0].conformer_id == "conf_001"
        assert conformers[0].energy_hartree == pytest.approx(-12.34567890)
        assert "3\n" in conformers[0].xyz_block
        assert "-12.34567890" in conformers[0].xyz_block
        assert "C  0.0  0.0  0.0" in conformers[0].xyz_block

    def test_multi_conformer(self, tmp_path):
        """parse_crest_ensemble correctly parses a multi-conformer file."""
        # Create test XYZ file with three conformers
        xyz_file = tmp_path / "multi.xyz"
        xyz_content = """3
-12.34567890
C  0.0  0.0  0.0
H  1.0  0.0  0.0
H  0.0  1.0  0.0
3
-12.34560000
C  0.1  0.0  0.0
H  1.1  0.0  0.0
H  0.1  1.0  0.0
3
-12.34555555
C  0.2  0.0  0.0
H  1.2  0.0  0.0
H  0.2  1.0  0.0
"""
        xyz_file.write_text(xyz_content)

        # Parse
        conformers = parse_crest_ensemble(xyz_file)

        # Verify count and IDs
        assert len(conformers) == 3
        assert conformers[0].conformer_id == "conf_001"
        assert conformers[1].conformer_id == "conf_002"
        assert conformers[2].conformer_id == "conf_003"

        # Verify energies
        assert conformers[0].energy_hartree == pytest.approx(-12.34567890)
        assert conformers[1].energy_hartree == pytest.approx(-12.34560000)
        assert conformers[2].energy_hartree == pytest.approx(-12.34555555)

        # Verify xyz_blocks are distinct
        assert "C  0.0  0.0  0.0" in conformers[0].xyz_block
        assert "C  0.1  0.0  0.0" in conformers[1].xyz_block
        assert "C  0.2  0.0  0.0" in conformers[2].xyz_block

    def test_xyz_block_structure(self, tmp_path):
        """parse_crest_ensemble reconstructs xyz_block with atom count, comment, and coordinates."""
        xyz_file = tmp_path / "test.xyz"
        xyz_content = """2
-5.5
O  0.0  0.0  0.0
H  1.0  0.0  0.0
"""
        xyz_file.write_text(xyz_content)

        conformers = parse_crest_ensemble(xyz_file)

        # xyz_block should contain all three components
        xyz_block = conformers[0].xyz_block
        lines = xyz_block.strip().split("\n")
        assert len(lines) == 4  # atom count + comment + 2 atoms
        assert lines[0].strip() == "2"
        assert "-5.5" in lines[1]
        assert "O  0.0  0.0  0.0" in lines[2]
        assert "H  1.0  0.0  0.0" in lines[3]


# Tests for CRESTConformer
class TestCRESTConformer:
    """Test CRESTConformer NamedTuple."""

    def test_structure(self):
        """CRESTConformer has correct fields."""
        conformer = CRESTConformer(
            conformer_id="conf_005",
            energy_hartree=-123.456,
            xyz_block="test content",
        )
        assert conformer.conformer_id == "conf_005"
        assert conformer.energy_hartree == -123.456
        assert conformer.xyz_block == "test content"


# Tests for ALPB_SOLVENT_MAP constant
class TestALPBSolventMap:
    """Test ALPB_SOLVENT_MAP constant."""

    def test_map_contents(self):
        """ALPB_SOLVENT_MAP contains expected solvent mappings."""
        assert ALPB_SOLVENT_MAP == {"chcl3": "chcl3", "dmso": "dmso"}
