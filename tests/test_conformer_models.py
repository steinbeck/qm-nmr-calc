"""Unit tests for conformer data models and storage."""

import tempfile
from pathlib import Path

import pytest

from qm_nmr_calc.models import (
    ConformerData,
    ConformerEnsemble,
    JobInput,
    JobStatus,
)
import qm_nmr_calc.storage as st
from qm_nmr_calc.storage import (
    create_conformer_directories,
    get_conformer_output_dir,
    get_conformer_scratch_dir,
    get_optimized_conformers_dir,
)


class TestConformerData:
    """Tests for ConformerData model."""

    def test_conformer_data_creation(self):
        """Create ConformerData with just conformer_id, check defaults."""
        c = ConformerData(conformer_id="conf_001")
        assert c.conformer_id == "conf_001"
        assert c.energy is None
        assert c.energy_unit is None
        assert c.weight is None
        assert c.geometry_file is None
        assert c.optimized_geometry_file is None
        assert c.rmsd_from_ref is None
        assert c.status == "pending"
        assert c.error_message is None

    def test_conformer_data_with_energy(self):
        """Create ConformerData with energy and energy_unit, verify values."""
        c = ConformerData(
            conformer_id="conf_002",
            energy=-154.321,
            energy_unit="hartree",
            weight=0.42,
        )
        assert c.conformer_id == "conf_002"
        assert c.energy == -154.321
        assert c.energy_unit == "hartree"
        assert c.weight == 0.42

    def test_conformer_data_status_transitions(self):
        """Verify all valid status values accepted."""
        valid_statuses = [
            "pending",
            "optimizing",
            "optimized",
            "nmr_running",
            "nmr_complete",
            "failed",
        ]
        for status in valid_statuses:
            c = ConformerData(conformer_id="conf_001", status=status)
            assert c.status == status

    def test_conformer_data_invalid_status(self):
        """Verify invalid status raises ValidationError."""
        from pydantic import ValidationError

        with pytest.raises(ValidationError):
            ConformerData(conformer_id="conf_001", status="invalid_status")

    def test_conformer_data_serialization(self):
        """Model_dump and model_validate roundtrip."""
        c1 = ConformerData(
            conformer_id="conf_003",
            energy=-200.5,
            energy_unit="kcal_mol",
            weight=0.15,
            geometry_file="output/conformers/conf_003.xyz",
            status="optimized",
        )
        dumped = c1.model_dump()
        c2 = ConformerData.model_validate(dumped)
        assert c2.conformer_id == c1.conformer_id
        assert c2.energy == c1.energy
        assert c2.energy_unit == c1.energy_unit
        assert c2.weight == c1.weight
        assert c2.geometry_file == c1.geometry_file
        assert c2.status == c1.status


class TestConformerEnsemble:
    """Tests for ConformerEnsemble model."""

    def test_ensemble_creation_minimal(self):
        """Method + empty conformers list."""
        e = ConformerEnsemble(method="rdkit_kdg", conformers=[])
        assert e.method == "rdkit_kdg"
        assert e.conformers == []
        assert e.temperature_k == 298.15
        assert e.pre_dft_energy_window_kcal == 6.0
        assert e.post_dft_energy_window_kcal == 3.0

    def test_ensemble_creation_with_conformers(self):
        """Method + list of ConformerData."""
        c1 = ConformerData(conformer_id="conf_001")
        c2 = ConformerData(conformer_id="conf_002")
        e = ConformerEnsemble(method="crest", conformers=[c1, c2])
        assert e.method == "crest"
        assert len(e.conformers) == 2
        assert e.conformers[0].conformer_id == "conf_001"
        assert e.conformers[1].conformer_id == "conf_002"

    def test_ensemble_defaults(self):
        """Temperature=298.15, pre_dft=6.0, post_dft=3.0."""
        e = ConformerEnsemble(method="rdkit_kdg", conformers=[])
        assert e.temperature_k == 298.15
        assert e.pre_dft_energy_window_kcal == 6.0
        assert e.post_dft_energy_window_kcal == 3.0
        assert e.total_generated == 0
        assert e.total_after_pre_filter == 0
        assert e.total_after_post_filter == 0

    def test_ensemble_custom_parameters(self):
        """Override temperature and energy windows."""
        e = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[],
            temperature_k=310.0,
            pre_dft_energy_window_kcal=10.0,
            post_dft_energy_window_kcal=5.0,
            total_generated=50,
            total_after_pre_filter=20,
            total_after_post_filter=8,
        )
        assert e.temperature_k == 310.0
        assert e.pre_dft_energy_window_kcal == 10.0
        assert e.post_dft_energy_window_kcal == 5.0
        assert e.total_generated == 50
        assert e.total_after_pre_filter == 20
        assert e.total_after_post_filter == 8

    def test_ensemble_serialization(self):
        """Model_dump/model_validate roundtrip with conformers."""
        c1 = ConformerData(conformer_id="conf_001", energy=-100.0, energy_unit="hartree")
        c2 = ConformerData(conformer_id="conf_002", energy=-99.8, energy_unit="hartree")
        e1 = ConformerEnsemble(
            method="rdkit_kdg",
            conformers=[c1, c2],
            temperature_k=298.15,
            total_generated=10,
        )
        dumped = e1.model_dump()
        e2 = ConformerEnsemble.model_validate(dumped)
        assert e2.method == e1.method
        assert len(e2.conformers) == len(e1.conformers)
        assert e2.conformers[0].conformer_id == "conf_001"
        assert e2.conformers[1].conformer_id == "conf_002"
        assert e2.temperature_k == e1.temperature_k
        assert e2.total_generated == e1.total_generated


class TestJobStatusBackwardCompatibility:
    """Tests for JobStatus backward compatibility with v1.x."""

    def test_jobstatus_defaults_single_conformer(self):
        """New JobStatus has conformer_mode='single', ensemble=None."""
        j = JobStatus(
            job_id="test123",
            status="queued",
            created_at="2024-01-01T00:00:00",
            input=JobInput(smiles="CCO", solvent="chcl3"),
            nwchem_version="7.0.2",
        )
        assert j.conformer_mode == "single"
        assert j.conformer_ensemble is None

    def test_jobstatus_with_ensemble(self):
        """Create JobStatus with conformer_ensemble populated."""
        c = ConformerData(conformer_id="conf_001")
        e = ConformerEnsemble(method="rdkit_kdg", conformers=[c])
        j = JobStatus(
            job_id="test456",
            status="running",
            created_at="2024-01-01T00:00:00",
            input=JobInput(smiles="CCCO", solvent="chcl3", conformer_mode="ensemble"),
            nwchem_version="7.0.2",
            conformer_mode="ensemble",
            conformer_ensemble=e,
        )
        assert j.conformer_mode == "ensemble"
        assert j.conformer_ensemble is not None
        assert j.conformer_ensemble.method == "rdkit_kdg"
        assert len(j.conformer_ensemble.conformers) == 1

    def test_jobstatus_loads_v1_json(self):
        """Create a dict mimicking v1.x status.json (no conformer fields), verify JobStatus.model_validate succeeds."""
        # Simulate v1.x status.json (no conformer_mode, no conformer_ensemble)
        v1_data = {
            "job_id": "abc123",
            "status": "complete",
            "created_at": "2024-01-01T00:00:00",
            "started_at": "2024-01-01T00:01:00",
            "completed_at": "2024-01-01T00:15:00",
            "input": {
                "smiles": "CCO",
                "solvent": "chcl3",
                "preset": "production",
            },
            "nwchem_version": "7.0.2",
            "nmr_results": None,
        }
        # This should load without error and apply defaults
        j = JobStatus.model_validate(v1_data)
        assert j.job_id == "abc123"
        assert j.status == "complete"
        assert j.conformer_mode == "single"  # Default
        assert j.conformer_ensemble is None  # Default


class TestJobInputConformerFields:
    """Tests for JobInput conformer fields."""

    def test_jobinput_defaults(self):
        """Conformer_mode='single', conformer_method=None."""
        ji = JobInput(smiles="CCO", solvent="chcl3")
        assert ji.conformer_mode == "single"
        assert ji.conformer_method is None
        assert ji.max_conformers is None

    def test_jobinput_ensemble_mode(self):
        """Conformer_mode='ensemble', conformer_method='rdkit_kdg'."""
        ji = JobInput(
            smiles="CCCO",
            solvent="dmso",
            conformer_mode="ensemble",
            conformer_method="rdkit_kdg",
            max_conformers=20,
        )
        assert ji.conformer_mode == "ensemble"
        assert ji.conformer_method == "rdkit_kdg"
        assert ji.max_conformers == 20


class TestStorageDirectories:
    """Tests for storage directory helpers."""

    @pytest.fixture
    def temp_data_dir(self, monkeypatch, tmp_path):
        """Monkeypatch DATA_DIR to use tmp_path."""
        temp_jobs = tmp_path / "jobs"
        temp_jobs.mkdir()
        monkeypatch.setattr(st, "DATA_DIR", temp_jobs)
        return temp_jobs

    def test_create_conformer_directories(self, temp_data_dir):
        """Creates correct directory structure."""
        job_id = "test_job_001"
        job_dir = temp_data_dir / job_id
        job_dir.mkdir()
        (job_dir / "scratch").mkdir()
        (job_dir / "output").mkdir()

        conformer_ids = ["conf_001", "conf_002", "conf_003"]
        result = create_conformer_directories(job_id, conformer_ids)

        assert len(result) == 3
        for cid in conformer_ids:
            assert cid in result
            scratch_dir = result[cid]
            assert scratch_dir.exists()
            assert scratch_dir == temp_data_dir / job_id / "scratch" / "conformers" / cid

            # Check output dir also created
            output_dir = temp_data_dir / job_id / "output" / "conformers" / cid
            assert output_dir.exists()

        # Check optimized dir created
        opt_dir = temp_data_dir / job_id / "output" / "optimized"
        assert opt_dir.exists()

    def test_conformer_scratch_isolation(self, temp_data_dir):
        """Each conformer gets unique scratch dir."""
        job_id = "test_job_002"
        job_dir = temp_data_dir / job_id
        job_dir.mkdir()
        (job_dir / "scratch").mkdir()
        (job_dir / "output").mkdir()

        conformer_ids = ["conf_001", "conf_002"]
        result = create_conformer_directories(job_id, conformer_ids)

        # Ensure scratch directories are isolated
        scratch1 = result["conf_001"]
        scratch2 = result["conf_002"]
        assert scratch1 != scratch2
        assert scratch1.parent == scratch2.parent  # Same parent
        assert scratch1.name == "conf_001"
        assert scratch2.name == "conf_002"

    def test_get_conformer_scratch_dir_path(self, temp_data_dir):
        """Returns correct path."""
        job_id = "test_job_003"
        job_dir = temp_data_dir / job_id
        job_dir.mkdir()
        (job_dir / "scratch").mkdir()
        (job_dir / "output").mkdir()

        create_conformer_directories(job_id, ["conf_001"])

        scratch_path = get_conformer_scratch_dir(job_id, "conf_001")
        assert scratch_path == temp_data_dir / job_id / "scratch" / "conformers" / "conf_001"
        assert scratch_path.exists()

    def test_get_conformer_output_dir_path(self, temp_data_dir):
        """Returns correct path."""
        job_id = "test_job_004"
        job_dir = temp_data_dir / job_id
        job_dir.mkdir()
        (job_dir / "scratch").mkdir()
        (job_dir / "output").mkdir()

        create_conformer_directories(job_id, ["conf_002"])

        output_path = get_conformer_output_dir(job_id, "conf_002")
        assert output_path == temp_data_dir / job_id / "output" / "conformers" / "conf_002"
        assert output_path.exists()

    def test_get_optimized_conformers_dir_path(self, temp_data_dir):
        """Returns correct path."""
        job_id = "test_job_005"
        job_dir = temp_data_dir / job_id
        job_dir.mkdir()
        (job_dir / "scratch").mkdir()
        (job_dir / "output").mkdir()

        create_conformer_directories(job_id, ["conf_001"])

        opt_path = get_optimized_conformers_dir(job_id)
        assert opt_path == temp_data_dir / job_id / "output" / "optimized"
        assert opt_path.exists()
