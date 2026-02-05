"""Job directory and status file management."""

import os
import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional

import orjson

from .models import JobInput, JobStatus, StepTiming

# Default data directory - relative to working directory
DATA_DIR = Path("./data/jobs")


def generate_job_id() -> str:
    """Generate URL-safe job ID (12 hex characters)."""
    return uuid.uuid4().hex[:12]


def get_job_dir(job_id: str) -> Path:
    """Get job directory path."""
    return DATA_DIR / job_id


def create_job_directory(
    smiles: str,
    solvent: str,
    nwchem_version: str,
    name: Optional[str] = None,
    preset: str = "production",
    notification_email: Optional[str] = None,
    conformer_mode: str = "single",
    conformer_method: Optional[str] = None,
    max_conformers: Optional[int] = None,
) -> JobStatus:
    """Create job directory with initial queued status.

    Creates a job directory structure:
        data/jobs/{job_id}/
            status.json  - Job metadata and status
            output/      - Calculation outputs
            logs/        - NWChem logs

    Args:
        smiles: SMILES string of molecule
        solvent: Solvent name for COSMO
        nwchem_version: NWChem version string
        name: Optional molecule name
        preset: Calculation preset (draft or production)
        notification_email: Optional email for completion notification
        conformer_mode: "single" (v1.x behavior) or "ensemble" (v2.0)
        conformer_method: Conformer generation method ("rdkit_kdg" or "crest"), only used when mode=ensemble
        max_conformers: Maximum conformers to generate (None = adaptive)

    Returns:
        JobStatus with initial queued status
    """
    job_id = generate_job_id()
    job_dir = get_job_dir(job_id)
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "output").mkdir(exist_ok=True)
    (job_dir / "logs").mkdir(exist_ok=True)

    # Make directories world-writable for multi-container setups
    # (API and worker may run as different UIDs, e.g., ARM64 local dev)
    os.chmod(job_dir, 0o777)
    os.chmod(job_dir / "output", 0o777)
    os.chmod(job_dir / "logs", 0o777)

    status = JobStatus(
        job_id=job_id,
        status="queued",
        created_at=datetime.utcnow(),
        input=JobInput(
            smiles=smiles,
            name=name,
            preset=preset,
            solvent=solvent,
            notification_email=notification_email,
            conformer_mode=conformer_mode,
            conformer_method=conformer_method,
            max_conformers=max_conformers,
        ),
        nwchem_version=nwchem_version,
        conformer_mode=conformer_mode,
    )

    _write_status(job_id, status)
    return status


def load_job_status(job_id: str) -> Optional[JobStatus]:
    """Load job status from disk.

    Returns None if job not found.
    """
    status_file = get_job_dir(job_id) / "status.json"
    if not status_file.exists():
        return None
    data = orjson.loads(status_file.read_bytes())
    return JobStatus.model_validate(data)


def update_job_status(job_id: str, **updates) -> JobStatus:
    """Update specific fields in job status.

    Raises ValueError if job not found.
    """
    status = load_job_status(job_id)
    if status is None:
        raise ValueError(f"Job {job_id} not found")
    updated = status.model_copy(update=updates)
    _write_status(job_id, updated)
    return updated


def start_step(job_id: str, step: str, label: str) -> JobStatus:
    """Start a new calculation step, completing the previous one if any.

    Args:
        job_id: Job identifier
        step: Step identifier (e.g., 'geometry_optimization')
        label: Human-readable label (e.g., 'Optimizing geometry')

    Returns:
        Updated JobStatus
    """
    status = load_job_status(job_id)
    if status is None:
        raise ValueError(f"Job {job_id} not found")

    now = datetime.utcnow()
    steps_completed = list(status.steps_completed)

    # Complete the previous step if one was in progress
    if status.current_step and status.step_started_at:
        duration = (now - status.step_started_at).total_seconds()
        steps_completed.append(
            StepTiming(
                step=status.current_step,
                label=status.current_step_label or status.current_step,
                started_at=status.step_started_at,
                completed_at=now,
                duration_seconds=round(duration, 1),
            )
        )

    updated = status.model_copy(
        update={
            "current_step": step,
            "current_step_label": label,
            "step_started_at": now,
            "steps_completed": steps_completed,
        }
    )
    _write_status(job_id, updated)
    return updated


def complete_current_step(job_id: str) -> JobStatus:
    """Complete the current step without starting a new one.

    Call this when the job is finishing.
    """
    status = load_job_status(job_id)
    if status is None:
        raise ValueError(f"Job {job_id} not found")

    now = datetime.utcnow()
    steps_completed = list(status.steps_completed)

    # Complete the current step if one was in progress
    if status.current_step and status.step_started_at:
        duration = (now - status.step_started_at).total_seconds()
        steps_completed.append(
            StepTiming(
                step=status.current_step,
                label=status.current_step_label or status.current_step,
                started_at=status.step_started_at,
                completed_at=now,
                duration_seconds=round(duration, 1),
            )
        )

    updated = status.model_copy(
        update={
            "current_step": None,
            "current_step_label": None,
            "step_started_at": None,
            "steps_completed": steps_completed,
        }
    )
    _write_status(job_id, updated)
    return updated


def _write_status(job_id: str, status: JobStatus) -> None:
    """Write status to disk as JSON."""
    status_file = get_job_dir(job_id) / "status.json"
    status_file.write_bytes(
        orjson.dumps(status.model_dump(mode="json"), option=orjson.OPT_INDENT_2)
    )
    # Make file world-writable for multi-container setups (different UIDs)
    # Only chmod if we're the owner (ignore errors when updating others' files)
    try:
        os.chmod(status_file, 0o666)
    except PermissionError:
        pass  # File already world-writable from original creation


def list_jobs_by_status(status_filter: str) -> list[str]:
    """List job IDs with given status.

    Returns empty list if no jobs match or data directory doesn't exist.
    """
    jobs = []
    if not DATA_DIR.exists():
        return jobs
    for job_dir in DATA_DIR.iterdir():
        if not job_dir.is_dir():
            continue
        status = load_job_status(job_dir.name)
        if status and status.status == status_filter:
            jobs.append(job_dir.name)
    return jobs


def get_geometry_file(job_id: str) -> Optional[Path]:
    """Get path to optimized geometry XYZ file if it exists."""
    geometry_file = get_job_dir(job_id) / "output" / "optimized.xyz"
    return geometry_file if geometry_file.exists() else None


def get_initial_geometry_file(job_id: str) -> Optional[Path]:
    """Get path to initial RDKit geometry XYZ file if it exists."""
    geometry_file = get_job_dir(job_id) / "output" / "initial.xyz"
    return geometry_file if geometry_file.exists() else None


def get_output_files(job_id: str) -> list[Path]:
    """Get list of raw NWChem output files in job scratch directory.

    Returns .out and .nw files from the scratch directory, including
    conformer subdirectories for ensemble calculations.
    """
    scratch_dir = get_job_dir(job_id) / "scratch"
    if not scratch_dir.exists():
        return []
    files = []
    # Top-level scratch files (single-conformer jobs)
    for pattern in ["*.out", "*.nw"]:
        files.extend(scratch_dir.glob(pattern))
    # Conformer subdirectory files (ensemble jobs)
    for pattern in ["*.out", "*.nw"]:
        files.extend(scratch_dir.glob(f"conformers/*/{pattern}"))
    return files


def get_visualization_file(job_id: str, filename: str) -> Optional[Path]:
    """Get path to a visualization file if it exists.

    Args:
        job_id: Job identifier
        filename: One of: spectrum_1H.svg, spectrum_1H.png, spectrum_13C.svg,
                  spectrum_13C.png, structure_annotated.svg, structure_annotated.png

    Returns:
        Path to file if exists, None otherwise
    """
    viz_file = get_job_dir(job_id) / "output" / filename
    return viz_file if viz_file.exists() else None


# v2.0 Conformer directory helpers


def create_conformer_directories(job_id: str, conformer_ids: list[str]) -> dict[str, Path]:
    """Create per-conformer scratch and output directories.

    Creates directory structure:
        data/jobs/{job_id}/
            scratch/conformers/{conf_id}/  - Isolated NWChem scratch per conformer
            output/conformers/{conf_id}/   - Per-conformer outputs
            output/optimized/              - Optimized geometries

    This prevents NWChem database file conflicts when running multiple conformers.

    Args:
        job_id: Job identifier
        conformer_ids: List of conformer IDs (e.g., ["conf_001", "conf_002"])

    Returns:
        Dict mapping conformer_id -> scratch_dir Path
    """
    job_dir = get_job_dir(job_id)
    scratch_base = job_dir / "scratch" / "conformers"
    output_base = job_dir / "output" / "conformers"
    optimized_dir = job_dir / "output" / "optimized"

    # Create parent directories
    scratch_base.mkdir(parents=True, exist_ok=True)
    output_base.mkdir(parents=True, exist_ok=True)
    optimized_dir.mkdir(parents=True, exist_ok=True)

    # Make parent directories world-writable for multi-container setups
    os.chmod(scratch_base, 0o777)
    os.chmod(output_base, 0o777)
    os.chmod(optimized_dir, 0o777)
    # Also chmod parents that were created
    os.chmod(job_dir / "scratch", 0o777)
    os.chmod(job_dir / "output" / "conformers", 0o777)

    # Create per-conformer directories
    result = {}
    for conf_id in conformer_ids:
        conf_scratch = scratch_base / conf_id
        conf_output = output_base / conf_id
        conf_scratch.mkdir(exist_ok=True)
        conf_output.mkdir(exist_ok=True)
        os.chmod(conf_scratch, 0o777)
        os.chmod(conf_output, 0o777)
        result[conf_id] = conf_scratch

    return result


def get_conformer_scratch_dir(job_id: str, conformer_id: str) -> Path:
    """Get path to conformer-specific scratch directory.

    Args:
        job_id: Job identifier
        conformer_id: Conformer ID (e.g., "conf_001")

    Returns:
        Path to scratch/conformers/{conformer_id}/

    Note: Does NOT create directory. Use create_conformer_directories first.
    """
    return get_job_dir(job_id) / "scratch" / "conformers" / conformer_id


def get_conformer_output_dir(job_id: str, conformer_id: str) -> Path:
    """Get path to conformer-specific output directory.

    Args:
        job_id: Job identifier
        conformer_id: Conformer ID (e.g., "conf_001")

    Returns:
        Path to output/conformers/{conformer_id}/
    """
    return get_job_dir(job_id) / "output" / "conformers" / conformer_id


def get_optimized_conformers_dir(job_id: str) -> Path:
    """Get path to optimized conformers directory.

    Args:
        job_id: Job identifier

    Returns:
        Path to output/optimized/
    """
    return get_job_dir(job_id) / "output" / "optimized"
