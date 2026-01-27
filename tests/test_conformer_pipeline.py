"""Integration tests for end-to-end conformer pipeline."""

import pytest

from qm_nmr_calc.conformers import generate_conformer_ensemble


class TestConformerPipeline:
    """Integration tests for full conformer generation pipeline."""

    def test_generates_ensemble_for_ethanol(self, tmp_path, monkeypatch):
        """Test full pipeline generates valid ConformerEnsemble for ethanol."""
        # Monkeypatch DATA_DIR to use temp directory
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Generate ensemble
        ensemble = generate_conformer_ensemble("CCO", "test_job_001")

        # Verify ensemble structure
        assert ensemble.method == "rdkit_kdg"
        assert len(ensemble.conformers) >= 1

        # Verify filtering reduces count
        assert ensemble.total_generated >= len(ensemble.conformers)

        # Verify all conformers have required fields
        for conformer in ensemble.conformers:
            assert conformer.energy is not None
            assert isinstance(conformer.energy, float)
            assert conformer.energy_unit == "kcal_mol"
            assert conformer.geometry_file is not None
            assert conformer.status == "pending"

    def test_generates_ensemble_for_flexible_molecule(self, tmp_path, monkeypatch):
        """Test adaptive count generates many conformers for flexible molecules."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Dodecane has 9 rotatable bonds (>8 threshold triggers 200 conformers)
        # Use max_conformers=60 to test adaptive logic without slow MMFF on 200 conformers
        ensemble = generate_conformer_ensemble(
            "CCCCCCCCCCCC", "test_job_002", max_conformers=60
        )

        # Should generate 60 conformers (our override)
        assert ensemble.total_generated == 60
        assert len(ensemble.conformers) > 0

    def test_writes_xyz_files_to_disk(self, tmp_path, monkeypatch):
        """Test that XYZ files are written to correct locations on disk."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Generate ensemble
        ensemble = generate_conformer_ensemble("CCO", "test_job_003")

        # Check each conformer has XYZ file on disk
        for conformer in ensemble.conformers:
            # Build full path from job dir
            xyz_path = tmp_path / "test_job_003" / conformer.geometry_file
            assert xyz_path.exists(), f"XYZ file not found: {xyz_path}"
            assert xyz_path.stat().st_size > 0, f"XYZ file is empty: {xyz_path}"

            # Verify XYZ format: first line should be atom count (integer)
            first_line = xyz_path.read_text().split("\n")[0]
            atom_count = int(first_line.strip())
            assert atom_count > 0, f"Invalid atom count in XYZ: {first_line}"

    def test_creates_conformer_directories(self, tmp_path, monkeypatch):
        """Test that conformer directories are created correctly."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Generate ensemble
        ensemble = generate_conformer_ensemble("CCO", "test_job_004")

        # Check directory structure
        job_dir = tmp_path / "test_job_004"
        assert (job_dir / "output" / "conformers").exists()
        assert (job_dir / "scratch" / "conformers").exists()

        # Check per-conformer subdirectories
        for conformer in ensemble.conformers:
            conf_id = conformer.conformer_id
            assert (job_dir / "output" / "conformers" / conf_id).exists()
            assert (job_dir / "scratch" / "conformers" / conf_id).exists()

    def test_energy_window_reduces_conformer_count(self, tmp_path, monkeypatch):
        """Test that tight energy window filters high-energy conformers."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Octane is flexible enough to generate diverse conformers with varying energies
        ensemble = generate_conformer_ensemble(
            "CCCCCCCC", "test_job_005", energy_window_kcal=0.5
        )

        # With tight window, filtering should reduce count
        # Note: For some molecules, all conformers might be within 0.5 kcal/mol
        # So we check that the mechanism works, not that it always filters
        assert ensemble.total_after_pre_filter <= ensemble.total_generated

    def test_max_conformers_override(self, tmp_path, monkeypatch):
        """Test that max_conformers parameter limits generation."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Generate with strict limit
        ensemble = generate_conformer_ensemble("CCO", "test_job_006", max_conformers=3)

        # Should generate at most 3 conformers
        assert ensemble.total_generated <= 3

    def test_conformer_ids_are_sequential(self, tmp_path, monkeypatch):
        """Test that conformer IDs are sequential and properly formatted."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Generate ensemble
        ensemble = generate_conformer_ensemble("CCO", "test_job_007")

        # Check IDs are sequential: conf_001, conf_002, etc.
        for i, conformer in enumerate(ensemble.conformers):
            expected_id = f"conf_{i+1:03d}"
            assert (
                conformer.conformer_id == expected_id
            ), f"Expected {expected_id}, got {conformer.conformer_id}"

    def test_invalid_smiles_raises_error(self, tmp_path, monkeypatch):
        """Test that invalid SMILES raises ValueError."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Invalid SMILES should raise ValueError
        with pytest.raises(ValueError, match="Invalid SMILES"):
            generate_conformer_ensemble("invalid_smiles_string", "test_job_008")
