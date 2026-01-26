"""Backward compatibility tests for v1.x JSON data with updated v2.0 models.

Verifies that existing v1.x status.json files (without conformer fields) load
correctly into the updated JobStatus/JobInput models, which now include optional
conformer-related fields with defaults.
"""

import pytest

from qm_nmr_calc.models import JobInput, JobStatus


# ---------------------------------------------------------------------------
# v1.x fixture data -- exact structure of a real v1.x status.json
# These dicts intentionally have NO conformer fields.
# ---------------------------------------------------------------------------

V1_STATUS_JSON = {
    "job_id": "abc123def456",
    "status": "complete",
    "created_at": "2026-01-20T10:00:00",
    "started_at": "2026-01-20T10:01:00",
    "completed_at": "2026-01-20T10:15:00",
    "input": {
        "smiles": "CCO",
        "name": "Ethanol",
        "preset": "production",
        "solvent": "chcl3",
        "notification_email": None,
    },
    "nwchem_version": "7.0.2",
    "current_step": None,
    "current_step_label": None,
    "step_started_at": None,
    "steps_completed": [
        {
            "step": "geometry_optimization",
            "label": "Optimizing geometry",
            "started_at": "2026-01-20T10:01:00",
            "completed_at": "2026-01-20T10:10:00",
            "duration_seconds": 540.0,
        }
    ],
    "error_message": None,
    "error_traceback": None,
    "cpu_time_seconds": 540.0,
    "memory_peak_mb": 256.0,
    "nmr_results": {
        "h1_shifts": [
            {"index": 7, "atom": "H", "shielding": 29.5, "shift": 3.7}
        ],
        "c13_shifts": [
            {"index": 1, "atom": "C", "shielding": 120.0, "shift": 58.0}
        ],
        "functional": "b3lyp",
        "basis_set": "6-311+G(2d,p)",
        "solvent": "chcl3",
    },
    "optimized_geometry_file": "output/optimized.xyz",
}

V1_QUEUED_STATUS = {
    "job_id": "queued999aaa",
    "status": "queued",
    "created_at": "2026-01-20T12:00:00",
    "input": {
        "smiles": "c1ccccc1",
        "preset": "draft",
        "solvent": "dmso",
    },
    "nwchem_version": "7.0.2",
}

V1_FAILED_STATUS = {
    "job_id": "fail456bbb00",
    "status": "failed",
    "created_at": "2026-01-20T14:00:00",
    "started_at": "2026-01-20T14:01:00",
    "completed_at": "2026-01-20T14:02:00",
    "input": {
        "smiles": "CC(=O)O",
        "name": "Acetic Acid",
        "preset": "production",
        "solvent": "chcl3",
    },
    "nwchem_version": "7.0.2",
    "error_message": "NWChem exited with code 1: SCF convergence failure",
    "error_traceback": "Traceback: ...",
}


class TestV1StatusLoading:
    """Verify v1.x status.json files load into updated JobStatus without error."""

    def test_v1_status_json_loads(self):
        """A v1.x status.json (no conformer fields) loads via model_validate."""
        job = JobStatus.model_validate(V1_STATUS_JSON)
        assert job.job_id == "abc123def456"
        assert job.status == "complete"

    def test_v1_status_defaults(self):
        """Loaded v1 job has conformer_mode='single' and conformer_ensemble=None."""
        job = JobStatus.model_validate(V1_STATUS_JSON)
        assert job.conformer_mode == "single"
        assert job.conformer_ensemble is None

    def test_v1_input_defaults(self):
        """Loaded v1 job's input has correct conformer defaults."""
        job = JobStatus.model_validate(V1_STATUS_JSON)
        assert job.input.conformer_mode == "single"
        assert job.input.conformer_method is None
        assert job.input.max_conformers is None

    def test_v1_status_preserves_all_fields(self):
        """Loaded v1 job has correct job_id, status, nmr_results -- no data corruption."""
        job = JobStatus.model_validate(V1_STATUS_JSON)

        # Core identification
        assert job.job_id == "abc123def456"
        assert job.status == "complete"
        assert job.nwchem_version == "7.0.2"

        # Input
        assert job.input.smiles == "CCO"
        assert job.input.name == "Ethanol"
        assert job.input.preset == "production"
        assert job.input.solvent == "chcl3"
        assert job.input.notification_email is None

        # Timestamps
        assert job.created_at.year == 2026
        assert job.started_at is not None
        assert job.completed_at is not None

        # Steps
        assert len(job.steps_completed) == 1
        assert job.steps_completed[0].step == "geometry_optimization"
        assert job.steps_completed[0].duration_seconds == 540.0

        # Resources
        assert job.cpu_time_seconds == 540.0
        assert job.memory_peak_mb == 256.0

        # NMR results
        assert job.nmr_results is not None
        assert len(job.nmr_results.h1_shifts) == 1
        assert job.nmr_results.h1_shifts[0].shift == 3.7
        assert len(job.nmr_results.c13_shifts) == 1
        assert job.nmr_results.c13_shifts[0].shift == 58.0
        assert job.nmr_results.functional == "b3lyp"
        assert job.nmr_results.basis_set == "6-311+G(2d,p)"
        assert job.nmr_results.solvent == "chcl3"

        # Geometry
        assert job.optimized_geometry_file == "output/optimized.xyz"

    def test_v1_status_roundtrip(self):
        """Load v1 JSON, dump to dict, reload -- all values preserved."""
        job1 = JobStatus.model_validate(V1_STATUS_JSON)
        dumped = job1.model_dump(mode="python")
        job2 = JobStatus.model_validate(dumped)

        assert job1.job_id == job2.job_id
        assert job1.status == job2.status
        assert job1.conformer_mode == job2.conformer_mode
        assert job1.conformer_ensemble == job2.conformer_ensemble
        assert job1.input.smiles == job2.input.smiles
        assert job1.input.conformer_mode == job2.input.conformer_mode
        assert job1.nmr_results == job2.nmr_results

    def test_v1_queued_job_loads(self):
        """Minimal v1 queued job (no results, no timestamps) loads correctly."""
        job = JobStatus.model_validate(V1_QUEUED_STATUS)
        assert job.job_id == "queued999aaa"
        assert job.status == "queued"
        assert job.started_at is None
        assert job.completed_at is None
        assert job.nmr_results is None
        assert job.conformer_mode == "single"
        assert job.conformer_ensemble is None

    def test_v1_failed_job_loads(self):
        """V1 failed job (with error_message) loads correctly."""
        job = JobStatus.model_validate(V1_FAILED_STATUS)
        assert job.job_id == "fail456bbb00"
        assert job.status == "failed"
        assert job.error_message == "NWChem exited with code 1: SCF convergence failure"
        assert job.error_traceback == "Traceback: ..."
        assert job.conformer_mode == "single"


class TestV1JobInputLoading:
    """Verify v1.x JobInput data loads correctly with conformer defaults."""

    def test_v1_job_input_loads(self):
        """Minimal v1.x JobInput loads with correct conformer defaults."""
        job_input = JobInput.model_validate(
            {"smiles": "CCO", "solvent": "chcl3"}
        )
        assert job_input.smiles == "CCO"
        assert job_input.solvent == "chcl3"
        assert job_input.conformer_mode == "single"
        assert job_input.conformer_method is None
        assert job_input.max_conformers is None
