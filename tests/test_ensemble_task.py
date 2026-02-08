"""Tests for ensemble NMR task and mode-based dispatch."""

import inspect
from unittest.mock import MagicMock, patch

import pytest


class TestEnsembleTaskImports:
    """Tests for ensemble task imports and basic structure."""

    def test_run_ensemble_nmr_task_imports(self):
        """Verify run_ensemble_nmr_task can be imported."""
        from qm_nmr_calc.tasks import run_ensemble_nmr_task
        assert callable(run_ensemble_nmr_task)

    def test_progress_callback_in_runner_signature(self):
        """Test that progress_callback is in run_ensemble_dft_and_nmr signature."""
        from qm_nmr_calc.nwchem.runner import run_ensemble_dft_and_nmr
        sig = inspect.signature(run_ensemble_dft_and_nmr)
        assert 'progress_callback' in sig.parameters

    def test_progress_callback_in_dft_optimization_signature(self):
        """Test that progress_callback is in run_conformer_dft_optimization signature."""
        from qm_nmr_calc.nwchem.runner import run_conformer_dft_optimization
        sig = inspect.signature(run_conformer_dft_optimization)
        assert 'progress_callback' in sig.parameters

    def test_progress_callback_in_nmr_calculations_signature(self):
        """Test that progress_callback is in run_conformer_nmr_calculations signature."""
        from qm_nmr_calc.nwchem.runner import run_conformer_nmr_calculations
        sig = inspect.signature(run_conformer_nmr_calculations)
        assert 'progress_callback' in sig.parameters


class TestStorageConformerParams:
    """Tests for storage.py conformer parameters."""

    def test_create_job_directory_accepts_conformer_mode(self):
        """Verify create_job_directory accepts conformer_mode parameter."""
        from qm_nmr_calc.storage import create_job_directory
        sig = inspect.signature(create_job_directory)
        assert 'conformer_mode' in sig.parameters

    def test_create_job_directory_accepts_conformer_method(self):
        """Verify create_job_directory accepts conformer_method parameter."""
        from qm_nmr_calc.storage import create_job_directory
        sig = inspect.signature(create_job_directory)
        assert 'conformer_method' in sig.parameters

    def test_create_job_directory_accepts_max_conformers(self):
        """Verify create_job_directory accepts max_conformers parameter."""
        from qm_nmr_calc.storage import create_job_directory
        sig = inspect.signature(create_job_directory)
        assert 'max_conformers' in sig.parameters


class TestEnsembleTaskFiltering:
    """Tests for ensemble task conformer filtering before averaging."""

    @patch('qm_nmr_calc.tasks.generate_conformer_ensemble')
    @patch('qm_nmr_calc.tasks.run_ensemble_dft_and_nmr')
    @patch('qm_nmr_calc.tasks.average_ensemble_nmr')
    @patch('qm_nmr_calc.tasks.load_job_status')
    @patch('qm_nmr_calc.tasks.update_job_status')
    @patch('qm_nmr_calc.tasks.start_step')
    @patch('qm_nmr_calc.tasks.complete_current_step')
    @patch('qm_nmr_calc.tasks.get_job_dir')
    @patch('qm_nmr_calc.tasks.generate_spectrum_plot')
    @patch('qm_nmr_calc.tasks.generate_annotated_structure')
    def test_ensemble_task_filters_nmr_complete_conformers(
        self, mock_annotate, mock_spectrum, mock_get_job_dir, mock_complete, mock_start, mock_update,
        mock_load, mock_avg, mock_dft, mock_gen
    ):
        """Verify ensemble task filters to nmr_complete conformers before averaging."""
        from qm_nmr_calc.tasks import run_ensemble_nmr_task
        from qm_nmr_calc.models import ConformerEnsemble, ConformerData, NMRResults, AtomShift

        # Setup mock job status
        mock_job = MagicMock()
        mock_job.status = 'queued'
        mock_job.input.smiles = 'CCO'
        mock_job.input.solvent = 'chcl3'
        mock_job.input.preset = 'production'
        mock_job.input.conformer_method = 'rdkit_kdg'
        mock_job.input.max_conformers = None
        mock_load.return_value = mock_job

        # Mock job directory
        mock_job_dir = MagicMock()
        mock_output_dir = MagicMock()
        mock_job_dir.__truediv__ = MagicMock(return_value=mock_output_dir)
        mock_get_job_dir.return_value = mock_job_dir

        # Setup ensemble with mixed statuses
        conformers = [
            ConformerData(conformer_id="conf_001", status="nmr_complete", energy=-100.0, energy_unit="hartree"),
            ConformerData(conformer_id="conf_002", status="failed"),
            ConformerData(conformer_id="conf_003", status="nmr_complete", energy=-99.5, energy_unit="hartree"),
        ]
        ensemble = ConformerEnsemble(method="rdkit_kdg", conformers=conformers)
        mock_gen.return_value = ensemble

        # NMR results only for nmr_complete (2 results)
        nmr1 = NMRResults(
            h1_shifts=[AtomShift(index=1, atom="H", shielding=30.0, shift=2.0)],
            c13_shifts=[AtomShift(index=2, atom="C", shielding=150.0, shift=50.0)],
            functional="b3lyp",
            basis_set="6-311+g(2d,p)",
            solvent="chcl3"
        )
        nmr2 = NMRResults(
            h1_shifts=[AtomShift(index=1, atom="H", shielding=31.0, shift=1.8)],
            c13_shifts=[AtomShift(index=2, atom="C", shielding=151.0, shift=49.0)],
            functional="b3lyp",
            basis_set="6-311+g(2d,p)",
            solvent="chcl3"
        )
        mock_dft.return_value = (ensemble, [nmr1, nmr2], None)

        # Mock avg_nmr return
        mock_avg.return_value = NMRResults(
            h1_shifts=[AtomShift(index=1, atom="H", shielding=30.5, shift=1.9)],
            c13_shifts=[AtomShift(index=2, atom="C", shielding=150.5, shift=49.5)],
            functional="b3lyp",
            basis_set="6-311+g(2d,p)",
            solvent="chcl3"
        )

        # Run task (calling the underlying function, not the Huey task wrapper)
        run_ensemble_nmr_task.call_local("test_job_id")

        # Verify average_ensemble_nmr was called with filtered ensemble (2 conformers, not 3)
        mock_avg.assert_called_once()
        filtered_ensemble = mock_avg.call_args[0][0]
        assert len(filtered_ensemble.conformers) == 2
        assert all(c.status == "nmr_complete" for c in filtered_ensemble.conformers)


class TestModeAwareDispatch:
    """Tests for mode-aware task dispatch in API endpoints."""

    @patch('qm_nmr_calc.api.routers.jobs.run_ensemble_nmr_task')
    @patch('qm_nmr_calc.api.routers.jobs.run_nmr_task')
    def test_submit_smiles_dispatches_ensemble_task(self, mock_single_task, mock_ensemble_task):
        """Verify submit_smiles dispatches to run_ensemble_nmr_task for ensemble mode."""
        from fastapi.testclient import TestClient
        from qm_nmr_calc.api.app import app

        client = TestClient(app)

        # Submit with ensemble mode
        response = client.post(
            "/api/v1/jobs",
            json={"smiles": "CCO", "solvent": "chcl3", "conformer_mode": "ensemble"},
        )

        # Should be accepted
        assert response.status_code == 202

        # Verify ensemble task was called, not single task
        mock_ensemble_task.assert_called_once()
        mock_single_task.assert_not_called()

    @patch('qm_nmr_calc.api.routers.jobs.run_ensemble_nmr_task')
    @patch('qm_nmr_calc.api.routers.jobs.run_nmr_task')
    def test_submit_smiles_dispatches_single_task(self, mock_single_task, mock_ensemble_task):
        """Verify submit_smiles dispatches to run_nmr_task for single mode."""
        from fastapi.testclient import TestClient
        from qm_nmr_calc.api.app import app

        client = TestClient(app)

        # Submit with single mode (explicitly)
        response = client.post(
            "/api/v1/jobs",
            json={"smiles": "CCO", "solvent": "chcl3", "conformer_mode": "single"},
        )

        # Should be accepted
        assert response.status_code == 202

        # Verify single task was called, not ensemble task
        mock_single_task.assert_called_once()
        mock_ensemble_task.assert_not_called()

    @patch('qm_nmr_calc.api.routers.jobs.run_ensemble_nmr_task')
    @patch('qm_nmr_calc.api.routers.jobs.run_nmr_task')
    def test_submit_smiles_returns_conformer_mode_in_response(self, mock_single_task, mock_ensemble_task):
        """Verify submit_smiles returns conformer_mode in response."""
        from fastapi.testclient import TestClient
        from qm_nmr_calc.api.app import app

        client = TestClient(app)

        # Submit with ensemble mode
        response = client.post(
            "/api/v1/jobs",
            json={"smiles": "CCO", "solvent": "chcl3", "conformer_mode": "ensemble"},
        )

        assert response.status_code == 202
        data = response.json()
        assert data["conformer_mode"] == "ensemble"


class TestConformerMethodWarning:
    """Tests for CREST fallback warning handling."""

    def test_job_status_has_conformer_method_warning_field(self):
        """Verify JobStatus model has conformer_method_warning field."""
        from qm_nmr_calc.models import JobStatus
        import inspect

        # Check that the field exists in model_fields
        assert 'conformer_method_warning' in JobStatus.model_fields

    @patch('qm_nmr_calc.tasks.generate_conformer_ensemble')
    @patch('qm_nmr_calc.tasks.run_ensemble_dft_and_nmr')
    @patch('qm_nmr_calc.tasks.average_ensemble_nmr')
    @patch('qm_nmr_calc.tasks.load_job_status')
    @patch('qm_nmr_calc.tasks.update_job_status')
    @patch('qm_nmr_calc.tasks.start_step')
    @patch('qm_nmr_calc.tasks.complete_current_step')
    @patch('qm_nmr_calc.tasks.get_job_dir')
    @patch('qm_nmr_calc.tasks.generate_spectrum_plot')
    @patch('qm_nmr_calc.tasks.generate_annotated_structure')
    def test_crest_fallback_sets_warning(
        self, mock_annotate, mock_spectrum, mock_get_job_dir, mock_complete, mock_start, mock_update,
        mock_load, mock_avg, mock_dft, mock_gen
    ):
        """Test CREST unavailable falls back to RDKit with warning."""
        from qm_nmr_calc.tasks import run_ensemble_nmr_task
        from qm_nmr_calc.models import ConformerEnsemble, ConformerData, NMRResults, AtomShift

        # Setup mock job status with crest method
        mock_job = MagicMock()
        mock_job.status = 'queued'
        mock_job.input.smiles = 'CCO'
        mock_job.input.solvent = 'chcl3'
        mock_job.input.preset = 'production'
        mock_job.input.conformer_method = 'crest'  # Request CREST
        mock_job.input.max_conformers = None
        mock_load.return_value = mock_job

        mock_job_dir = MagicMock()
        mock_output_dir = MagicMock()
        mock_job_dir.__truediv__ = MagicMock(return_value=mock_output_dir)
        mock_get_job_dir.return_value = mock_job_dir

        # First call raises CREST unavailable error, second succeeds with RDKit
        conformers = [
            ConformerData(conformer_id="conf_001", status="nmr_complete", energy=-100.0, energy_unit="hartree"),
        ]
        ensemble = ConformerEnsemble(method="rdkit_kdg", conformers=conformers)

        mock_gen.side_effect = [
            ValueError("CREST conformer method requested but CREST/xTB not installed."),
            ensemble
        ]

        nmr1 = NMRResults(
            h1_shifts=[AtomShift(index=1, atom="H", shielding=30.0, shift=2.0)],
            c13_shifts=[AtomShift(index=2, atom="C", shielding=150.0, shift=50.0)],
            functional="b3lyp",
            basis_set="6-311+g(2d,p)",
            solvent="chcl3"
        )
        mock_dft.return_value = (ensemble, [nmr1], None)
        mock_avg.return_value = nmr1

        # Run task
        run_ensemble_nmr_task.call_local("test_job_id")

        # Verify update_job_status was called with warning
        warning_calls = [
            call for call in mock_update.call_args_list
            if 'conformer_method_warning' in call.kwargs
        ]
        assert len(warning_calls) >= 1
        warning_msg = warning_calls[0].kwargs['conformer_method_warning']
        assert "CREST" in warning_msg
        assert "RDKit" in warning_msg
