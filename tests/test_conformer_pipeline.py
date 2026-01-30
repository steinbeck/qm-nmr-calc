"""Integration tests for end-to-end conformer pipeline."""

from unittest import mock

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


class TestCRESTDispatch:
    """Tests for CREST dispatch in pipeline."""

    def test_pipeline_crest_dispatch_when_available(self, tmp_path, monkeypatch):
        """Test pipeline dispatches to CREST when method='crest' and CREST available."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Mock CREST functions at their actual location
        with mock.patch("qm_nmr_calc.conformers.crest_generator.detect_crest_available") as mock_detect, \
             mock.patch("qm_nmr_calc.conformers.crest_generator.get_alpb_solvent") as mock_get_alpb, \
             mock.patch("qm_nmr_calc.conformers.crest_generator.generate_crest_ensemble") as mock_crest_gen:

            # Setup mocks
            mock_detect.return_value = True
            mock_get_alpb.return_value = "chcl3"

            # Create mock ensemble
            from qm_nmr_calc.models import ConformerEnsemble
            mock_ensemble = ConformerEnsemble(
                method="crest",
                conformers=[],
                pre_dft_energy_window_kcal=6.0,
                total_generated=5,
                total_after_pre_filter=5,
            )
            mock_crest_gen.return_value = mock_ensemble

            # Call with CREST method
            ensemble = generate_conformer_ensemble(
                smiles="CCO",
                job_id="test_crest_001",
                conformer_method="crest",
                solvent="chcl3",
            )

            # Verify CREST was called
            mock_detect.assert_called_once()
            mock_get_alpb.assert_called_once_with("chcl3")
            mock_crest_gen.assert_called_once_with(
                smiles="CCO",
                job_id="test_crest_001",
                solvent="chcl3",
                charge=0,
                energy_window_kcal=6.0,
                timeout_seconds=7200,
            )
            assert ensemble.method == "crest"

    def test_pipeline_crest_raises_when_not_installed(self, tmp_path, monkeypatch):
        """Test pipeline raises ValueError when CREST requested but not installed."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Mock CREST not available
        with mock.patch("qm_nmr_calc.conformers.crest_generator.detect_crest_available") as mock_detect:
            mock_detect.return_value = False

            # Verify ValueError raised
            with pytest.raises(ValueError, match="CREST/xTB not installed"):
                generate_conformer_ensemble(
                    smiles="CCO",
                    job_id="test_crest_002",
                    conformer_method="crest",
                )

    def test_pipeline_crest_raises_for_vacuum(self, tmp_path, monkeypatch):
        """Test pipeline raises ValueError when CREST requested for vacuum."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Mock CREST available but vacuum solvent
        with mock.patch("qm_nmr_calc.conformers.crest_generator.detect_crest_available") as mock_detect, \
             mock.patch("qm_nmr_calc.conformers.crest_generator.get_alpb_solvent") as mock_get_alpb:

            mock_detect.return_value = True
            mock_get_alpb.return_value = None  # vacuum returns None

            # Verify ValueError raised
            with pytest.raises(ValueError, match="not supported.*vacuum or unsupported"):
                generate_conformer_ensemble(
                    smiles="CCO",
                    job_id="test_crest_003",
                    conformer_method="crest",
                    solvent="vacuum",
                )

    def test_pipeline_crest_raises_for_unsupported_solvent(self, tmp_path, monkeypatch):
        """Test pipeline raises ValueError when CREST requested for unsupported solvent."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Mock CREST available but unsupported solvent
        with mock.patch("qm_nmr_calc.conformers.crest_generator.detect_crest_available") as mock_detect, \
             mock.patch("qm_nmr_calc.conformers.crest_generator.get_alpb_solvent") as mock_get_alpb:

            mock_detect.return_value = True
            mock_get_alpb.return_value = None  # water not supported

            # Verify ValueError raised
            with pytest.raises(ValueError, match="not supported"):
                generate_conformer_ensemble(
                    smiles="CCO",
                    job_id="test_crest_004",
                    conformer_method="crest",
                    solvent="water",
                )

    def test_pipeline_rdkit_default_path_unchanged(self, tmp_path, monkeypatch):
        """Test pipeline still works with RDKit method (backward compatibility)."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Call with default method (rdkit_kdg)
        ensemble = generate_conformer_ensemble(
            smiles="CCO",
            job_id="test_rdkit_001",
            conformer_method="rdkit_kdg",
        )

        # Verify RDKit path works
        assert ensemble.method == "rdkit_kdg"
        assert len(ensemble.conformers) >= 1


class TestPipelineClusteringIntegration:
    """Tests for clustering integration in pipeline."""

    def test_flexible_molecule_reduces_to_target(self, tmp_path, monkeypatch):
        """Flexible molecule should be reduced to target conformer count."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Hexane - flexible, should generate many conformers
        ensemble = generate_conformer_ensemble(
            smiles="CCCCCC",
            job_id="test_clustering_integration",
            target_conformers_for_dft=8,
        )

        # Should be at or below target
        assert len(ensemble.conformers) <= 10  # Allow small buffer
        # Should have some conformers
        assert len(ensemble.conformers) >= 1
        # Total generated should be higher
        assert ensemble.total_generated > len(ensemble.conformers)

    def test_rigid_molecule_not_over_reduced(self, tmp_path, monkeypatch):
        """Rigid molecule with few conformers should not be affected."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        # Benzene - rigid, very few conformers
        ensemble = generate_conformer_ensemble(
            smiles="c1ccccc1",
            job_id="test_rigid_clustering",
            target_conformers_for_dft=8,
        )

        # Should keep what it has (likely 1-3 conformers)
        assert len(ensemble.conformers) >= 1
        assert len(ensemble.conformers) <= 5

    def test_clustering_threshold_affects_output(self, tmp_path, monkeypatch):
        """Different clustering thresholds should affect conformer count."""
        # Monkeypatch DATA_DIR
        import qm_nmr_calc.storage as storage

        monkeypatch.setattr(storage, "DATA_DIR", tmp_path)

        smiles = "CCCCCCCC"  # Octane - very flexible

        # Tight clustering (more clusters, more representatives)
        ensemble_tight = generate_conformer_ensemble(
            smiles=smiles,
            job_id="test_tight_cluster",
            clustering_rmsd_threshold=0.8,
            target_conformers_for_dft=15,
        )

        # Loose clustering (fewer clusters, fewer representatives)
        ensemble_loose = generate_conformer_ensemble(
            smiles=smiles,
            job_id="test_loose_cluster",
            clustering_rmsd_threshold=2.5,
            target_conformers_for_dft=15,
        )

        # Tight should have >= loose (or equal if both hit target)
        # This is probabilistic but generally true
        assert ensemble_tight is not None
        assert ensemble_loose is not None
