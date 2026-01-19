"""Job directory and status file management."""

import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional

import orjson

from .models import JobInput, JobStatus

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
    isicle_version: str,
    nwchem_version: str,
    name: Optional[str] = None,
    preset: str = "production",
) -> JobStatus:
    """Create job directory with initial queued status.

    Creates a job directory structure:
        data/jobs/{job_id}/
            status.json  - Job metadata and status
            output/      - Calculation outputs
            logs/        - NWChem logs
    """
    job_id = generate_job_id()
    job_dir = get_job_dir(job_id)
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "output").mkdir(exist_ok=True)
    (job_dir / "logs").mkdir(exist_ok=True)

    status = JobStatus(
        job_id=job_id,
        status="queued",
        created_at=datetime.utcnow(),
        input=JobInput(smiles=smiles, name=name, preset=preset, solvent=solvent),
        isicle_version=isicle_version,
        nwchem_version=nwchem_version,
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


def _write_status(job_id: str, status: JobStatus) -> None:
    """Write status to disk as JSON."""
    status_file = get_job_dir(job_id) / "status.json"
    status_file.write_bytes(
        orjson.dumps(status.model_dump(mode="json"), option=orjson.OPT_INDENT_2)
    )


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


def get_output_files(job_id: str) -> list[Path]:
    """Get list of raw NWChem output files in job scratch directory.

    Returns .out and .nw files from the scratch directory.
    """
    scratch_dir = get_job_dir(job_id) / "scratch"
    if not scratch_dir.exists():
        return []
    files = []
    for pattern in ["*.out", "*.nw"]:
        files.extend(scratch_dir.glob(pattern))
    return files
