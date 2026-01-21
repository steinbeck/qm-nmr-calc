"""Huey task definitions for calculations."""
from pathlib import Path

import orjson

from .queue import huey
from .storage import load_job_status, get_job_dir, update_job_status, start_step, complete_current_step
from .nwchem import run_calculation
from .presets import PRESETS, PresetName
from .shifts import shielding_to_shift
from .models import NMRResults, AtomShift
from .visualization import generate_spectrum_plot, generate_annotated_structure


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

    # Accept both 'queued' and 'running' - signal handler sets 'running' before task body
    if job_status.status not in ('queued', 'running'):
        raise ValueError(f"Job {job_id} is not ready to run (status: {job_status.status})")

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


@huey.task()
def run_nmr_task(job_id: str) -> dict:
    """
    Execute NMR calculation (geometry optimization + shielding) for a queued job.

    This is the main task for NMR calculations. It:
    1. Loads job status and validates it's queued
    2. Gets preset configuration based on job input
    3. Runs two-step DFT: geometry optimization then NMR shielding (both with COSMO)
    4. Converts shielding to chemical shifts
    5. Saves NMR results to job output directory
    6. Updates job status with NMR results

    Step progress is tracked and updated throughout.

    Args:
        job_id: ID of the job to process

    Returns:
        dict with success status and output file paths

    Note: Status updates (running/complete/failed) happen via signal handlers
    in queue.py. Let exceptions propagate for SIGNAL_ERROR.
    """
    # Load job info
    job_status = load_job_status(job_id)
    if job_status is None:
        raise ValueError(f"Job {job_id} not found")

    # Accept both 'queued' and 'running' - signal handler sets 'running' before task body
    if job_status.status not in ('queued', 'running'):
        raise ValueError(f"Job {job_id} is not ready to run (status: {job_status.status})")

    # Get preset configuration
    preset_name = PresetName(job_status.input.preset)
    preset = PRESETS[preset_name]

    smiles = job_status.input.smiles
    solvent = job_status.input.solvent
    job_dir = get_job_dir(job_id)

    # Step 1: Geometry optimization
    start_step(job_id, "geometry_optimization", "Optimizing geometry")

    # Step 2: NMR shielding calculation (tracked within run_calculation)
    start_step(job_id, "nmr_shielding", "Computing NMR shielding")

    # Run complete calculation (geometry optimization + NMR shielding)
    # Both steps use COSMO solvation - this fixes the gas-phase bug
    result = run_calculation(
        smiles=smiles,
        job_dir=job_dir,
        preset=preset,
        solvent=solvent,
        processes=preset['processes'],
    )

    geometry_file = result['geometry_file']

    # Step 3: Post-processing
    start_step(job_id, "post_processing", "Generating results")

    # Convert shielding to chemical shifts using preset-specific TMS reference
    shifts = shielding_to_shift(result['shielding_data'], preset=job_status.input.preset)

    # Build NMRResults object
    h1_shifts = [
        AtomShift(
            index=s['index'],
            atom=s['atom'],
            shielding=s['shielding'],
            shift=s['shift']
        )
        for s in shifts['1H']
    ]
    c13_shifts = [
        AtomShift(
            index=s['index'],
            atom=s['atom'],
            shielding=s['shielding'],
            shift=s['shift']
        )
        for s in shifts['13C']
    ]

    nmr_results = NMRResults(
        h1_shifts=h1_shifts,
        c13_shifts=c13_shifts,
        functional=preset['functional'],
        basis_set=preset['nmr_basis_set'],
        solvent=solvent,
    )

    # Save nmr_results.json to output directory
    output_dir = job_dir / 'output'
    results_file = output_dir / 'nmr_results.json'
    results_file.write_bytes(
        orjson.dumps(
            nmr_results.model_dump(mode='json'),
            option=orjson.OPT_INDENT_2
        )
    )

    # Generate visualizations
    # 1H spectrum
    generate_spectrum_plot(
        shifts=[s.shift for s in h1_shifts],
        nucleus="1H",
        output_dir=output_dir,
    )

    # 13C spectrum
    generate_spectrum_plot(
        shifts=[s.shift for s in c13_shifts],
        nucleus="13C",
        output_dir=output_dir,
    )

    # Annotated structure (both 1H and 13C on one image)
    generate_annotated_structure(
        smiles=smiles,
        h1_shifts=[{"index": s.index, "shift": s.shift} for s in h1_shifts],
        c13_shifts=[{"index": s.index, "shift": s.shift} for s in c13_shifts],
        output_dir=output_dir,
    )

    # Complete final step
    complete_current_step(job_id)

    # Update job status with NMR results
    update_job_status(
        job_id,
        nmr_results=nmr_results,
        optimized_geometry_file=str(geometry_file),
    )

    return {
        'success': True,
        'job_id': job_id,
        'output_files': {
            'geometry': str(geometry_file),
            'nmr_results': str(results_file),
        }
    }
