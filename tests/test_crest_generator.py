"""Tests for CREST conformer generation utilities."""

import os
import subprocess
import unittest.mock as mock
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from qm_nmr_calc.conformers.crest_generator import (
    ALPB_SOLVENT_MAP,
    CRESTConformer,
    detect_crest_available,
    generate_crest_ensemble,
    get_alpb_solvent,
    parse_crest_ensemble,
    run_crest,
)
from qm_nmr_calc.models import ConformerEnsemble


# Tests for detect_crest_available()
class TestDetectCRESTAvailable:
    """Test CREST/xTB binary detection."""

    def test_both_binaries_present(self):
        """detect_crest_available returns True when both crest and xtb are found."""
        # Clear cache before test
        detect_crest_available.cache_clear()
        with mock.patch("qm_nmr_calc.conformers.crest_generator.shutil.which") as mock_which:
            # Both binaries found
            mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}" if cmd in ["crest", "xtb"] else None
            result = detect_crest_available()
            assert result is True

    def test_only_crest_present(self):
        """detect_crest_available returns False when only crest is found."""
        # Clear cache before test
        detect_crest_available.cache_clear()
        with mock.patch("qm_nmr_calc.conformers.crest_generator.shutil.which") as mock_which:
            # Only crest found
            mock_which.side_effect = lambda cmd: "/usr/bin/crest" if cmd == "crest" else None
            result = detect_crest_available()
            assert result is False

    def test_only_xtb_present(self):
        """detect_crest_available returns False when only xtb is found."""
        # Clear cache before test
        detect_crest_available.cache_clear()
        with mock.patch("qm_nmr_calc.conformers.crest_generator.shutil.which") as mock_which:
            # Only xtb found
            mock_which.side_effect = lambda cmd: "/usr/bin/xtb" if cmd == "xtb" else None
            result = detect_crest_available()
            assert result is False

    def test_neither_binary_present(self):
        """detect_crest_available returns False when neither binary is found."""
        # Clear cache before test
        detect_crest_available.cache_clear()
        with mock.patch("qm_nmr_calc.conformers.crest_generator.shutil.which") as mock_which:
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


# Tests for run_crest()
class TestRunCrest:
    """Test CREST subprocess execution."""

    def test_success_returns_output_path(self, tmp_path):
        """Successful CREST run returns path to crest_conformers.xyz."""
        input_xyz = tmp_path / "crest_input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        output_file = tmp_path / "crest_conformers.xyz"
        output_file.write_text("3\n-12.345\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        with patch("qm_nmr_calc.conformers.crest_generator.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            result = run_crest(input_xyz, "chcl3", charge=0, ewin_kcal=6.0, timeout_seconds=60)

            assert result == output_file
            assert result.exists()
            mock_run.assert_called_once()

    def test_timeout_raises_runtime_error(self, tmp_path):
        """CREST timeout raises RuntimeError with helpful message."""
        input_xyz = tmp_path / "crest_input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        with patch("qm_nmr_calc.conformers.crest_generator.subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.TimeoutExpired(
                cmd=["crest"], timeout=60
            )

            with pytest.raises(RuntimeError) as exc_info:
                run_crest(input_xyz, "chcl3", timeout_seconds=60)

            assert "timeout" in str(exc_info.value).lower()
            assert "rdkit" in str(exc_info.value).lower()

    def test_nonzero_exit_raises_runtime_error(self, tmp_path):
        """CREST non-zero exit code raises RuntimeError."""
        input_xyz = tmp_path / "crest_input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        with patch("qm_nmr_calc.conformers.crest_generator.subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(
                returncode=1,
                cmd=["crest"],
                stderr="CREST error message"
            )

            with pytest.raises(RuntimeError) as exc_info:
                run_crest(input_xyz, "chcl3")

            assert "exit code 1" in str(exc_info.value)
            assert "CREST error message" in str(exc_info.value)

    def test_missing_output_file_raises_runtime_error(self, tmp_path):
        """CREST success but missing output file raises RuntimeError."""
        input_xyz = tmp_path / "crest_input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        # Don't create output file

        with patch("qm_nmr_calc.conformers.crest_generator.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            with pytest.raises(RuntimeError) as exc_info:
                run_crest(input_xyz, "chcl3")

            assert "crest_conformers.xyz" in str(exc_info.value)

    def test_command_args_built_correctly(self, tmp_path):
        """Verify CREST command arguments are correct."""
        input_xyz = tmp_path / "crest_input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        output_file = tmp_path / "crest_conformers.xyz"
        output_file.write_text("3\n-12.345\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        with patch("qm_nmr_calc.conformers.crest_generator.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            run_crest(
                input_xyz,
                "dmso",
                charge=-1,
                ewin_kcal=8.0,
                num_threads=4,
                timeout_seconds=120
            )

            call_args = mock_run.call_args
            cmd = call_args[0][0]

            assert cmd[0] == "crest"
            assert str(input_xyz) in cmd
            assert "--gfn2" in cmd
            assert "--alpb" in cmd
            assert "dmso" in cmd
            assert "--chrg" in cmd
            assert "-1" in cmd
            assert "--ewin" in cmd
            assert "8.0" in cmd
            assert "-T" in cmd
            assert "4" in cmd

            # Check keyword args
            assert call_args[1]["cwd"] == tmp_path
            assert call_args[1]["timeout"] == 120

    def test_environment_variables_set(self, tmp_path):
        """Verify OMP_STACKSIZE and GFORTRAN_UNBUFFERED_ALL are set."""
        input_xyz = tmp_path / "crest_input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        output_file = tmp_path / "crest_conformers.xyz"
        output_file.write_text("3\n-12.345\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        with patch("qm_nmr_calc.conformers.crest_generator.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            run_crest(input_xyz, "chcl3")

            call_args = mock_run.call_args
            env = call_args[1]["env"]

            assert env["OMP_STACKSIZE"] == "2G"
            assert env["GFORTRAN_UNBUFFERED_ALL"] == "1"

    def test_default_num_threads(self, tmp_path):
        """Default num_threads is os.cpu_count() or 4."""
        input_xyz = tmp_path / "crest_input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        output_file = tmp_path / "crest_conformers.xyz"
        output_file.write_text("3\n-12.345\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        with patch("qm_nmr_calc.conformers.crest_generator.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            run_crest(input_xyz, "chcl3")

            call_args = mock_run.call_args
            cmd = call_args[0][0]

            expected_threads = str(os.cpu_count() or 4)
            assert expected_threads in cmd


# Tests for generate_crest_ensemble()
class TestGenerateCrestEnsemble:
    """Test CREST ensemble generation pipeline."""

    def test_successful_ensemble_generation(self, tmp_path):
        """Successful ensemble generation produces correct ConformerEnsemble."""
        job_id = "test_job"

        # Mock run_crest to return fake ensemble file
        fake_ensemble = tmp_path / "crest_conformers.xyz"
        fake_ensemble.write_text(
            "3\n-12.345\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
            "3\n-12.343\nC 0.1 0.0 0.0\nH 1.1 0.0 0.0\nH 0.1 1.0 0.0\n"
        )

        with (
            patch("qm_nmr_calc.conformers.crest_generator.run_crest") as mock_run_crest,
            patch("qm_nmr_calc.conformers.crest_generator.get_job_dir") as mock_job_dir,
            patch("qm_nmr_calc.conformers.crest_generator.create_conformer_directories") as mock_create_dirs,
        ):
            mock_run_crest.return_value = fake_ensemble
            mock_job_dir.return_value = tmp_path
            mock_create_dirs.return_value = {}

            ensemble = generate_crest_ensemble(
                smiles="C",
                job_id=job_id,
                solvent="chcl3",
                charge=0,
                energy_window_kcal=6.0,
                timeout_seconds=60
            )

            assert isinstance(ensemble, ConformerEnsemble)
            assert ensemble.method == "crest"
            assert len(ensemble.conformers) == 2
            assert ensemble.total_generated == 2
            assert ensemble.total_after_pre_filter == 2

    def test_energy_conversion_hartree_to_kcal(self, tmp_path):
        """Energy conversion from Hartree to relative kcal/mol is correct."""
        job_id = "test_job"

        # Create ensemble with known energies
        # Min energy: -12.345 Hartree
        # Second energy: -12.343 Hartree (0.002 Hartree higher)
        # Expected relative: 0.0 and 1.255 kcal/mol (0.002 * 627.5095)
        fake_ensemble = tmp_path / "crest_conformers.xyz"
        fake_ensemble.write_text(
            "3\n-12.345\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
            "3\n-12.343\nC 0.1 0.0 0.0\nH 1.1 0.0 0.0\nH 0.1 1.0 0.0\n"
        )

        with (
            patch("qm_nmr_calc.conformers.crest_generator.run_crest") as mock_run_crest,
            patch("qm_nmr_calc.conformers.crest_generator.get_job_dir") as mock_job_dir,
            patch("qm_nmr_calc.conformers.crest_generator.create_conformer_directories") as mock_create_dirs,
        ):
            mock_run_crest.return_value = fake_ensemble
            mock_job_dir.return_value = tmp_path
            mock_create_dirs.return_value = {}

            ensemble = generate_crest_ensemble(
                smiles="C",
                job_id=job_id,
                solvent="chcl3"
            )

            # Check relative energies
            assert ensemble.conformers[0].energy == pytest.approx(0.0, abs=0.001)
            assert ensemble.conformers[1].energy == pytest.approx(1.255, abs=0.001)
            assert all(c.energy_unit == "kcal_mol" for c in ensemble.conformers)

    def test_energy_window_filter_applied(self, tmp_path):
        """Energy window filter removes high-energy conformers."""
        job_id = "test_job"

        # Create ensemble with one conformer outside energy window
        # Energy difference: 10 kcal/mol = 0.01594 Hartree
        fake_ensemble = tmp_path / "crest_conformers.xyz"
        fake_ensemble.write_text(
            "3\n-12.345\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
            "3\n-12.343\nC 0.1 0.0 0.0\nH 1.1 0.0 0.0\nH 0.1 1.0 0.0\n"
            "3\n-12.329\nC 0.2 0.0 0.0\nH 1.2 0.0 0.0\nH 0.2 1.0 0.0\n"
        )

        with (
            patch("qm_nmr_calc.conformers.crest_generator.run_crest") as mock_run_crest,
            patch("qm_nmr_calc.conformers.crest_generator.get_job_dir") as mock_job_dir,
            patch("qm_nmr_calc.conformers.crest_generator.create_conformer_directories") as mock_create_dirs,
        ):
            mock_run_crest.return_value = fake_ensemble
            mock_job_dir.return_value = tmp_path
            mock_create_dirs.return_value = {}

            ensemble = generate_crest_ensemble(
                smiles="C",
                job_id=job_id,
                solvent="chcl3",
                energy_window_kcal=6.0
            )

            # Third conformer should be filtered out (10 kcal/mol > 6.0 window)
            assert len(ensemble.conformers) == 2
            assert ensemble.total_generated == 3
            assert ensemble.total_after_pre_filter == 2

    def test_xyz_files_written_correctly(self, tmp_path):
        """Individual XYZ files written to correct paths."""
        job_id = "test_job"

        fake_ensemble = tmp_path / "crest_conformers.xyz"
        fake_ensemble.write_text(
            "3\n-12.345\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
            "3\n-12.343\nC 0.1 0.0 0.0\nH 1.1 0.0 0.0\nH 0.1 1.0 0.0\n"
        )

        with (
            patch("qm_nmr_calc.conformers.crest_generator.run_crest") as mock_run_crest,
            patch("qm_nmr_calc.conformers.crest_generator.get_job_dir") as mock_job_dir,
            patch("qm_nmr_calc.conformers.crest_generator.create_conformer_directories") as mock_create_dirs,
            patch("qm_nmr_calc.conformers.crest_generator.get_conformer_output_dir") as mock_output_dir,
        ):
            mock_run_crest.return_value = fake_ensemble
            mock_job_dir.return_value = tmp_path

            # Create actual directories
            conf_001_dir = tmp_path / "output" / "conformers" / "conf_001"
            conf_002_dir = tmp_path / "output" / "conformers" / "conf_002"
            conf_001_dir.mkdir(parents=True)
            conf_002_dir.mkdir(parents=True)

            def get_output_dir_side_effect(job_id, conf_id):
                return tmp_path / "output" / "conformers" / conf_id

            mock_output_dir.side_effect = get_output_dir_side_effect
            mock_create_dirs.return_value = {}

            ensemble = generate_crest_ensemble(
                smiles="C",
                job_id=job_id,
                solvent="chcl3"
            )

            # Verify XYZ files written
            xyz_001 = conf_001_dir / "geometry.xyz"
            xyz_002 = conf_002_dir / "geometry.xyz"

            assert xyz_001.exists()
            assert xyz_002.exists()

            # Verify content
            content_001 = xyz_001.read_text()
            assert "3\n" in content_001
            assert "-12.345" in content_001
            assert "C 0.0 0.0 0.0" in content_001

    def test_conformer_data_fields_populated(self, tmp_path):
        """ConformerData fields populated correctly."""
        job_id = "test_job"

        fake_ensemble = tmp_path / "crest_conformers.xyz"
        fake_ensemble.write_text(
            "3\n-12.345\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
            "3\n-12.343\nC 0.1 0.0 0.0\nH 1.1 0.0 0.0\nH 0.1 1.0 0.0\n"
        )

        with (
            patch("qm_nmr_calc.conformers.crest_generator.run_crest") as mock_run_crest,
            patch("qm_nmr_calc.conformers.crest_generator.get_job_dir") as mock_job_dir,
            patch("qm_nmr_calc.conformers.crest_generator.create_conformer_directories") as mock_create_dirs,
            patch("qm_nmr_calc.conformers.crest_generator.get_conformer_output_dir") as mock_output_dir,
        ):
            mock_run_crest.return_value = fake_ensemble
            mock_job_dir.return_value = tmp_path
            mock_create_dirs.return_value = {}

            def get_output_dir_side_effect(job_id, conf_id):
                out_dir = tmp_path / "output" / "conformers" / conf_id
                out_dir.mkdir(parents=True, exist_ok=True)
                return out_dir

            mock_output_dir.side_effect = get_output_dir_side_effect

            ensemble = generate_crest_ensemble(
                smiles="C",
                job_id=job_id,
                solvent="chcl3"
            )

            # Check first conformer
            conf_001 = ensemble.conformers[0]
            assert conf_001.conformer_id == "conf_001"
            assert conf_001.energy == pytest.approx(0.0, abs=0.001)
            assert conf_001.energy_unit == "kcal_mol"
            assert conf_001.geometry_file == "output/conformers/conf_001/geometry.xyz"
            assert conf_001.status == "pending"

            # Check second conformer
            conf_002 = ensemble.conformers[1]
            assert conf_002.conformer_id == "conf_002"
            assert conf_002.energy == pytest.approx(1.255, abs=0.001)
            assert conf_002.energy_unit == "kcal_mol"
            assert conf_002.geometry_file == "output/conformers/conf_002/geometry.xyz"
            assert conf_002.status == "pending"
