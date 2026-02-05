"""Huey task definitions for calculations."""
import os
from pathlib import Path

import orjson
from rdkit import Chem
from rdkit.Chem import AllChem

from .queue import huey
from .storage import load_job_status, get_job_dir, update_job_status, start_step, complete_current_step
from .nwchem import run_calculation
from .nwchem.runner import run_ensemble_dft_and_nmr
from .conformers.pipeline import generate_conformer_ensemble
from .conformers.boltzmann import average_ensemble_nmr
from .presets import PRESETS, PresetName
from .shifts import shielding_to_shift
from .models import NMRResults, AtomShift
from .visualization import generate_spectrum_plot, generate_annotated_structure


def _generate_initial_xyz(smiles: str, output_path: Path) -> None:
    """Generate and save RDKit 3D geometry for immediate visualization."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # ETKDGv3 with deterministic seed (project standard: 0xF00D)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMoleculeConfs(mol)

    # Generate XYZ
    xyz_lines = [str(mol.GetNumAtoms()), smiles]
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz_lines.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

    output_path.write_text("\n".join(xyz_lines))
    # Make file world-writable for multi-container setups (different UIDs)
    try:
        os.chmod(output_path, 0o666)
    except PermissionError:
        pass  # File created by another container with correct permissions


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

    # Generate initial RDKit geometry for 3D visualization (before optimization starts)
    # Skip if API already created it (multi-container permission handling)
    initial_xyz_path = job_dir / "output" / "initial.xyz"
    if not initial_xyz_path.exists():
        _generate_initial_xyz(smiles, initial_xyz_path)

    # Step 1: Geometry optimization
    start_step(job_id, "geometry_optimization", "Optimizing geometry")

    # Callback to switch step tracking when optimization completes
    def on_opt_complete():
        start_step(job_id, "nmr_shielding", "Computing NMR shielding")

    # Run complete calculation (geometry optimization + NMR shielding)
    # Both steps use COSMO solvation - this fixes the gas-phase bug
    result = run_calculation(
        smiles=smiles,
        job_dir=job_dir,
        preset=preset,
        solvent=solvent,
        processes=preset['processes'],
        on_optimization_complete=on_opt_complete,
    )

    geometry_file = result['geometry_file']

    # Store NWChem runtime info (actual processes and memory used)
    runtime_info = result.get('runtime_info')
    if runtime_info:
        update_job_status(
            job_id,
            nwchem_actual_processes=runtime_info['nproc'],
            nwchem_actual_memory_mb=runtime_info['total_memory_mb'],
        )

    # Step 3: Post-processing
    start_step(job_id, "post_processing", "Generating results")

    # Convert shielding to chemical shifts using DELTA50 regression factors
    shifts = shielding_to_shift(
        result['shielding_data'],
        functional=preset['functional'].upper(),  # 'b3lyp' -> 'B3LYP'
        basis_set=preset['nmr_basis_set'],
        solvent=solvent,
    )

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


@huey.task()
def run_ensemble_nmr_task(job_id: str) -> dict:
    """
    Execute ensemble NMR calculation (conformer generation + DFT + NMR + Boltzmann averaging).

    This is the v2.0 task for multi-conformer NMR calculations. It:
    1. Loads job status and validates it's queued
    2. Gets preset configuration based on job input
    3. Generates conformer ensemble (RDKit KDG or CREST)
    4. Runs DFT optimization on all conformers
    5. Applies post-DFT energy filter
    6. Runs NMR shielding on filtered conformers
    7. Computes Boltzmann-weighted average shifts
    8. Generates visualizations
    9. Updates job status with averaged results

    Args:
        job_id: ID of the job to process

    Returns:
        dict with success status

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

    # Determine conformer method (default to rdkit_kdg)
    conformer_method = job_status.input.conformer_method or "rdkit_kdg"

    # Step 1: Generate conformer ensemble
    start_step(job_id, "generating_conformers", "Generating conformers")
    try:
        ensemble = generate_conformer_ensemble(
            smiles=smiles,
            job_id=job_id,
            conformer_method=conformer_method,
            solvent=solvent,
            max_conformers=job_status.input.max_conformers,
        )
    except ValueError as e:
        # CREST unavailable or unsupported solvent - fall back to RDKit with warning
        if "CREST" in str(e) and conformer_method == "crest":
            warning_msg = str(e) + " Falling back to RDKit KDG method."
            update_job_status(job_id, conformer_method_warning=warning_msg)
            ensemble = generate_conformer_ensemble(
                smiles=smiles,
                job_id=job_id,
                conformer_method="rdkit_kdg",
                solvent=solvent,
                max_conformers=job_status.input.max_conformers,
            )
        else:
            raise

    update_job_status(job_id, conformer_ensemble=ensemble)

    # Progress callback to update job status with X/N counts
    def on_progress(step: str, current: int, total: int):
        """Update job status with progress during conformer processing."""
        start_step(job_id, step, f"{step.replace('_', ' ').title()} ({current}/{total})")

    # Step 2-4: DFT optimization + Post-DFT filter + NMR calculations
    start_step(job_id, "optimizing_conformers", f"Optimizing conformers (0/{len(ensemble.conformers)})")
    ensemble, nmr_results, runtime_info = run_ensemble_dft_and_nmr(
        ensemble=ensemble,
        job_id=job_id,
        preset=preset,
        solvent=solvent,
        processes=preset['processes'],
        progress_callback=on_progress,
    )

    # Store NWChem runtime info (actual processes and memory used)
    if runtime_info:
        update_job_status(
            job_id,
            nwchem_actual_processes=runtime_info['nproc'],
            nwchem_actual_memory_mb=runtime_info['total_memory_mb'],
        )

    # Step 5: Boltzmann averaging
    # CRITICAL: Filter to only nmr_complete conformers before averaging
    start_step(job_id, "averaging_shifts", "Computing Boltzmann average")
    nmr_complete_conformers = [c for c in ensemble.conformers if c.status == "nmr_complete"]

    if not nmr_complete_conformers:
        raise RuntimeError("No conformers completed NMR calculation successfully")

    # Create filtered ensemble for averaging (same count as nmr_results)
    filtered_ensemble = ensemble.model_copy(update={"conformers": nmr_complete_conformers})
    avg_nmr = average_ensemble_nmr(filtered_ensemble, nmr_results)

    # Step 6: Post-processing (visualizations)
    start_step(job_id, "post_processing", "Generating results")
    output_dir = job_dir / 'output'

    # Save nmr_results.json to output directory
    results_file = output_dir / 'nmr_results.json'
    results_file.write_bytes(
        orjson.dumps(
            avg_nmr.model_dump(mode='json'),
            option=orjson.OPT_INDENT_2
        )
    )

    # Generate spectrum plots using averaged shifts
    generate_spectrum_plot(
        shifts=[s.shift for s in avg_nmr.h1_shifts],
        nucleus="1H",
        output_dir=output_dir,
    )
    generate_spectrum_plot(
        shifts=[s.shift for s in avg_nmr.c13_shifts],
        nucleus="13C",
        output_dir=output_dir,
    )

    # Generate annotated structure using averaged shifts
    generate_annotated_structure(
        smiles=smiles,
        h1_shifts=[{"index": s.index, "shift": s.shift} for s in avg_nmr.h1_shifts],
        c13_shifts=[{"index": s.index, "shift": s.shift} for s in avg_nmr.c13_shifts],
        output_dir=output_dir,
    )

    complete_current_step(job_id)

    # Update job status with averaged NMR results and full ensemble data
    update_job_status(
        job_id,
        nmr_results=avg_nmr,
        conformer_ensemble=ensemble,  # Full ensemble with all conformers and their statuses
    )

    return {
        'success': True,
        'job_id': job_id,
    }
