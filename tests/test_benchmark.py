"""Tests for DELTA50 benchmark infrastructure."""

from collections import deque
from datetime import datetime, timezone
from pathlib import Path

import orjson
import pytest

from qm_nmr_calc.benchmark import (
    FAILURE_THRESHOLD,
    BenchmarkResult,
    ExperimentalShifts,
    MoleculeData,
    build_task_matrix,
    check_stop_requested,
    clear_stop_file,
    get_data_dir,
    load_delta50_molecules,
    load_experimental_shifts,
    update_status,
)
from qm_nmr_calc.benchmark.runner import (
    BENCHMARK_PRESETS,
    FUNCTIONALS,
    SOLVENTS,
    is_task_complete,
)


class TestDataLoader:
    """Tests for data_loader module."""

    def test_get_data_dir_exists(self):
        """Data directory should exist."""
        data_dir = get_data_dir()
        assert data_dir.exists(), f"Data directory not found: {data_dir}"

    def test_load_experimental_shifts(self):
        """Should load all 50 molecules with experimental data."""
        shifts = load_experimental_shifts()
        assert isinstance(shifts, ExperimentalShifts)
        assert (
            len(shifts.molecules) == 50
        ), f"Expected 50 molecules, got {len(shifts.molecules)}"
        assert "doi" in shifts.source

    def test_load_experimental_shifts_has_h1_data(self):
        """Each molecule should have H1 shift data."""
        shifts = load_experimental_shifts()
        for mol_id, mol in shifts.molecules.items():
            assert len(mol.h1_shifts) > 0, f"Molecule {mol_id} has no 1H shifts"

    def test_load_experimental_shifts_has_c13_data(self):
        """Each molecule should have C13 shift data."""
        shifts = load_experimental_shifts()
        for mol_id, mol in shifts.molecules.items():
            assert len(mol.c13_shifts) > 0, f"Molecule {mol_id} has no 13C shifts"

    def test_load_delta50_molecules(self):
        """Should load all 50 molecule structures."""
        molecules = load_delta50_molecules()
        assert len(molecules) == 50, f"Expected 50 molecules, got {len(molecules)}"

    def test_load_delta50_molecules_returns_rdkit_mol(self):
        """Each molecule should be a valid RDKit mol with XYZ block."""
        molecules = load_delta50_molecules()
        for mol_id, (mol, xyz_block) in molecules.items():
            assert mol is not None, f"Molecule {mol_id} is None"
            assert mol.GetNumAtoms() > 0, f"Molecule {mol_id} has no atoms"
            assert isinstance(xyz_block, str), f"XYZ block for {mol_id} is not string"


class TestModels:
    """Tests for Pydantic models."""

    def test_molecule_data_creation(self):
        """MoleculeData should validate correctly."""
        mol = MoleculeData(
            id="test_01",
            name="Test Molecule",
            xyz_file="molecules/test_01.xyz",
            h1_shifts=[7.26, 3.45],
            c13_shifts=[128.5, 45.2],
        )
        assert mol.id == "test_01"
        assert len(mol.h1_shifts) == 2

    def test_benchmark_result_creation(self):
        """BenchmarkResult should validate correctly."""
        result = BenchmarkResult(
            molecule_id="test_01",
            functional="B3LYP",
            basis_set="6-311+G(2d,p)",
            solvent="CHCl3",
            calculated_h1=[7.30, 3.50],
            calculated_c13=[130.0, 46.0],
            status="complete",
        )
        assert result.molecule_id == "test_01"
        assert result.status == "complete"


class TestRunner:
    """Tests for benchmark runner."""

    def test_benchmark_presets_defined(self):
        """Both B3LYP and WP04 presets should be defined."""
        assert "B3LYP" in BENCHMARK_PRESETS
        assert "WP04" in BENCHMARK_PRESETS

    def test_benchmark_presets_have_required_keys(self):
        """Presets should have all required configuration keys."""
        required_keys = {"functional", "basis_set", "nmr_basis_set", "max_iter"}
        for name, preset in BENCHMARK_PRESETS.items():
            missing = required_keys - set(preset.keys())
            assert not missing, f"Preset {name} missing keys: {missing}"

    def test_build_task_matrix_default(self):
        """Default matrix should be 50 molecules x 2 functionals x 2 solvents."""
        tasks = build_task_matrix()
        expected = 50 * len(FUNCTIONALS) * len(SOLVENTS)
        assert len(tasks) == expected, f"Expected {expected} tasks, got {len(tasks)}"

    def test_build_task_matrix_filtered(self):
        """Should filter by molecules, functionals, solvents."""
        tasks = build_task_matrix(
            molecules=["compound_01", "compound_02"],
            functionals=["B3LYP"],
            solvents=["CHCl3"],
        )
        assert len(tasks) == 2, f"Expected 2 tasks, got {len(tasks)}"

    def test_task_has_required_keys(self):
        """Each task should have task_id, molecule_id, functional, solvent."""
        tasks = build_task_matrix()
        for task in tasks[:5]:  # Check first 5
            assert "task_id" in task
            assert "molecule_id" in task
            assert "functional" in task
            assert "solvent" in task

    def test_is_task_complete_false_when_no_dir(self, tmp_path):
        """Task should not be complete if output dir doesn't exist."""
        task = {"molecule_id": "test", "functional": "B3LYP", "solvent": "CHCl3"}
        assert not is_task_complete(tmp_path, task)

    def test_is_task_complete_false_when_no_shifts_json(self, tmp_path):
        """Task should not be complete if shifts.json doesn't exist."""
        task = {"molecule_id": "test", "functional": "B3LYP", "solvent": "CHCl3"}
        output_dir = tmp_path / "test" / "B3LYP_CHCl3"
        output_dir.mkdir(parents=True)
        assert not is_task_complete(tmp_path, task)

    def test_is_task_complete_true_when_shifts_json_exists(self, tmp_path):
        """Task should be complete if shifts.json exists."""
        task = {"molecule_id": "test", "functional": "B3LYP", "solvent": "CHCl3"}
        output_dir = tmp_path / "test" / "B3LYP_CHCl3"
        output_dir.mkdir(parents=True)
        (output_dir / "shifts.json").write_text("{}")
        assert is_task_complete(tmp_path, task)


class TestStatusTracking:
    """Tests for status file tracking and control mechanisms."""

    def test_failure_threshold_defined(self):
        """Failure threshold should be 5 (10% of 50 molecules)."""
        assert FAILURE_THRESHOLD == 5

    def test_check_stop_requested_false_when_no_file(self, tmp_path):
        """Should return False when STOP file doesn't exist."""
        assert not check_stop_requested(tmp_path)

    def test_check_stop_requested_true_when_file_exists(self, tmp_path):
        """Should return True when STOP file exists."""
        (tmp_path / "STOP").touch()
        assert check_stop_requested(tmp_path)

    def test_clear_stop_file_removes_file(self, tmp_path):
        """Should remove STOP file if it exists."""
        stop_file = tmp_path / "STOP"
        stop_file.touch()
        assert stop_file.exists()
        clear_stop_file(tmp_path)
        assert not stop_file.exists()

    def test_clear_stop_file_no_error_if_missing(self, tmp_path):
        """Should not raise error if STOP file doesn't exist."""
        clear_stop_file(tmp_path)  # Should not raise

    def test_update_status_creates_file(self, tmp_path):
        """update_status should create status.json file."""
        started = datetime.now(timezone.utc)
        update_status(
            results_dir=tmp_path,
            state="running",
            total_tasks=200,
            completed=10,
            failed=2,
            current_task=None,
            failures=[],
            started_at=started,
            calc_times=deque([]),
        )
        status_file = tmp_path / "status.json"
        assert status_file.exists()

    def test_update_status_has_required_fields(self, tmp_path):
        """Status file should contain all required fields."""
        started = datetime.now(timezone.utc)
        update_status(
            results_dir=tmp_path,
            state="running",
            total_tasks=200,
            completed=10,
            failed=2,
            current_task={"molecule_id": "compound_01", "functional": "B3LYP", "solvent": "CHCl3"},
            failures=[{"molecule_id": "compound_02", "error": "test error"}],
            started_at=started,
            calc_times=deque([120.5, 130.2, 140.8]),
        )
        status_file = tmp_path / "status.json"
        data = orjson.loads(status_file.read_bytes())

        # Required fields
        assert "run_id" in data
        assert "started_at" in data
        assert "updated_at" in data
        assert data["state"] == "running"
        assert data["total_tasks"] == 200
        assert data["completed"] == 10
        assert data["failed"] == 2
        assert data["current_task"] is not None
        assert "failures" in data
        assert "estimated_remaining_hours" in data
        assert "avg_calc_time_seconds" in data

    def test_update_status_calculates_eta_from_rolling_average(self, tmp_path):
        """ETA should be calculated from rolling average of recent calc times."""
        started = datetime.now(timezone.utc)
        # 3 calculations averaging 120 seconds
        calc_times = deque([100.0, 120.0, 140.0])
        update_status(
            results_dir=tmp_path,
            state="running",
            total_tasks=200,
            completed=3,
            failed=0,
            current_task=None,
            failures=[],
            started_at=started,
            calc_times=calc_times,
        )
        status_file = tmp_path / "status.json"
        data = orjson.loads(status_file.read_bytes())

        # Average is 120s, 197 remaining, so 197 * 120 / 3600 = 6.57 hours
        assert data["avg_calc_time_seconds"] == 120.0
        assert data["estimated_remaining_hours"] is not None
        assert data["estimated_remaining_hours"] > 6.0

    def test_update_status_keeps_last_10_failures(self, tmp_path):
        """Should only keep last 10 failures."""
        started = datetime.now(timezone.utc)
        failures = [{"molecule_id": f"mol_{i}", "error": "test"} for i in range(15)]
        update_status(
            results_dir=tmp_path,
            state="running",
            total_tasks=200,
            completed=10,
            failed=15,
            current_task=None,
            failures=failures,
            started_at=started,
            calc_times=deque([]),
        )
        status_file = tmp_path / "status.json"
        data = orjson.loads(status_file.read_bytes())

        assert len(data["failures"]) == 10
        # Should be last 10 (mol_5 through mol_14)
        assert data["failures"][0]["molecule_id"] == "mol_5"

    def test_update_status_atomic_write(self, tmp_path):
        """Status file should not have .tmp suffix after write."""
        started = datetime.now(timezone.utc)
        update_status(
            results_dir=tmp_path,
            state="running",
            total_tasks=200,
            completed=0,
            failed=0,
            current_task=None,
            failures=[],
            started_at=started,
            calc_times=deque([]),
        )
        # No temp file should remain
        temp_file = tmp_path / "status.json.tmp"
        assert not temp_file.exists()
        # Final file should exist and be valid JSON
        status_file = tmp_path / "status.json"
        assert status_file.exists()
        orjson.loads(status_file.read_bytes())  # Should not raise


class TestExpandedSolvents:
    """Tests for expanded solvent list in benchmark runner."""

    def test_solvents_list_has_six_entries(self):
        """SOLVENTS list should have exactly 6 entries."""
        assert len(SOLVENTS) == 6, f"Expected 6 solvents, got {len(SOLVENTS)}"
        expected = {"CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene"}
        assert set(SOLVENTS) == expected

    def test_solvents_include_new_entries(self):
        """All 4 new solvents should be in SOLVENTS list."""
        assert "Methanol" in SOLVENTS
        assert "Water" in SOLVENTS
        assert "Acetone" in SOLVENTS
        assert "Benzene" in SOLVENTS

    def test_solvents_retain_original_entries(self):
        """Original CHCl3 and DMSO should still be present."""
        assert "CHCl3" in SOLVENTS
        assert "DMSO" in SOLVENTS

    def test_build_task_matrix_with_new_solvent(self):
        """Task matrix with single new solvent should produce correct tasks."""
        tasks = build_task_matrix(
            molecules=["compound_01"],
            functionals=["B3LYP"],
            solvents=["Benzene"],
        )
        assert len(tasks) == 1
        assert tasks[0]["solvent"] == "Benzene"
        assert tasks[0]["molecule_id"] == "compound_01"
        assert tasks[0]["functional"] == "B3LYP"

    def test_build_task_matrix_with_all_new_solvents(self):
        """Task matrix with all new solvents should produce correct count."""
        tasks = build_task_matrix(
            molecules=["compound_01"],
            functionals=["B3LYP"],
            solvents=["Methanol", "Water", "Acetone", "Benzene"],
        )
        assert len(tasks) == 4
        solvents_in_tasks = [t["solvent"] for t in tasks]
        assert "Methanol" in solvents_in_tasks
        assert "Water" in solvents_in_tasks
        assert "Acetone" in solvents_in_tasks
        assert "Benzene" in solvents_in_tasks

    def test_build_task_matrix_default_includes_new_solvents(self):
        """Default task matrix should include all 6 solvents."""
        tasks = build_task_matrix(
            molecules=["compound_01"],
            functionals=["B3LYP"],
        )
        assert len(tasks) == 6
        solvents_in_tasks = {t["solvent"] for t in tasks}
        assert solvents_in_tasks == {"CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene"}
