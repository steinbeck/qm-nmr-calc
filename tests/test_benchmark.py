"""Tests for DELTA50 benchmark infrastructure."""

import json
from pathlib import Path

import pytest

from qm_nmr_calc.benchmark import (
    BenchmarkResult,
    ExperimentalShifts,
    MoleculeData,
    build_task_matrix,
    get_data_dir,
    load_delta50_molecules,
    load_experimental_shifts,
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
