"""Huey task definitions for calculations."""
from pathlib import Path

from .queue import huey
from .storage import load_job_status, get_job_dir
from .isicle_wrapper import run_geometry_optimization


@huey.task()
def run_optimization_task(job_id: str) -> dict:
    """
    Execute geometry optimization for a queued job.

    Args:
        job_id: ID of the job to process

    Returns:
        dict with success status and output path

    Note: Status updates happen via signal handlers in queue.py,
    not in this function. Let exceptions propagate for SIGNAL_ERROR.
    """
    # Load job info
    job_status = load_job_status(job_id)
    if job_status is None:
        raise ValueError(f"Job {job_id} not found")

    if job_status.status != 'queued':
        raise ValueError(f"Job {job_id} is not queued (status: {job_status.status})")

    smiles = job_status.input.smiles
    job_dir = get_job_dir(job_id)

    # Run the calculation
    # Exceptions propagate to Huey -> SIGNAL_ERROR -> status update
    output_file = run_geometry_optimization(
        smiles=smiles,
        job_dir=job_dir,
        processes=4,  # Single worker, so can use multiple MPI processes
    )

    return {
        'success': True,
        'job_id': job_id,
        'output_file': str(output_file)
    }
