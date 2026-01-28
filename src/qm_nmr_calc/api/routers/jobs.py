"""Job submission and status endpoints."""

import io
import zipfile
from typing import Annotated, Optional

from fastapi import APIRouter, Form, HTTPException, UploadFile, status
from fastapi.responses import FileResponse, JSONResponse, Response
from rdkit import Chem

from ..schemas import JobStatusResponse, JobSubmitRequest, NMRResultsResponse, ProblemDetail
from ...nwchem import get_nwchem_version
from ...models import JobStatus
from ...shifts import get_scaling_factor
from ...solvents import validate_solvent, get_supported_solvents
from ...storage import create_job_directory, get_job_dir, get_geometry_file, get_initial_geometry_file, get_output_files, get_visualization_file, load_job_status
from ...tasks import run_nmr_task, run_ensemble_nmr_task, _generate_initial_xyz
from ...validation import validate_mol_file, validate_smiles

router = APIRouter(prefix="/jobs", tags=["jobs"])


def job_status_to_response(job_status: JobStatus) -> dict:
    """Convert JobStatus model to JobStatusResponse dict.

    Flattens input.smiles -> input_smiles, input.name -> input_name.
    Converts datetime fields to ISO strings.
    Converts NMR results if present.
    Builds ensemble-specific fields for conformer_mode="ensemble".
    """
    from ...conformers.boltzmann import HARTREE_TO_KCAL

    # Convert NMR results if present
    nmr_results = None
    if job_status.nmr_results is not None:
        # Get factor metadata for MAE display
        # Note: functional stored lowercase, get_scaling_factor expects uppercase
        h1_factor = get_scaling_factor(
            job_status.nmr_results.functional.upper(),
            job_status.nmr_results.basis_set,
            "1H",
            job_status.nmr_results.solvent,
        )
        c13_factor = get_scaling_factor(
            job_status.nmr_results.functional.upper(),
            job_status.nmr_results.basis_set,
            "13C",
            job_status.nmr_results.solvent,
        )

        nmr_results = {
            "h1_shifts": [
                {"index": s.index, "atom": s.atom, "shift": s.shift}
                for s in job_status.nmr_results.h1_shifts
            ],
            "c13_shifts": [
                {"index": s.index, "atom": s.atom, "shift": s.shift}
                for s in job_status.nmr_results.c13_shifts
            ],
            "functional": job_status.nmr_results.functional,
            "basis_set": job_status.nmr_results.basis_set,
            "solvent": job_status.nmr_results.solvent,
            "scaling_factor_source": "DELTA50",
            "h1_expected_mae": f"+/- {h1_factor['mae']:.2f} ppm",
            "c13_expected_mae": f"+/- {c13_factor['mae']:.2f} ppm",
        }

    # Convert completed steps to response format
    steps_completed = [
        {
            "step": step.step,
            "label": step.label,
            "started_at": step.started_at.isoformat() + "Z",
            "completed_at": step.completed_at.isoformat() + "Z",
            "duration_seconds": step.duration_seconds,
        }
        for step in job_status.steps_completed
    ]

    # Build ensemble-specific fields
    conformer_method = None
    conformer_count = None
    conformer_progress = None
    ensemble_metadata = None
    conformer_method_warning = getattr(job_status, "conformer_method_warning", None)

    if job_status.conformer_mode == "ensemble" and job_status.conformer_ensemble:
        ensemble = job_status.conformer_ensemble
        conformer_method = ensemble.method
        conformer_count = len(ensemble.conformers)

        # Find minimum energy for relative energy display
        min_energy = None
        for c in ensemble.conformers:
            if c.energy is not None and c.status in ("optimized", "nmr_complete"):
                if min_energy is None or c.energy < min_energy:
                    min_energy = c.energy

        # Build per-conformer progress array
        conformer_progress = []
        for c in ensemble.conformers:
            energy_kcal = None
            if c.energy is not None and min_energy is not None:
                if c.energy_unit == "hartree":
                    energy_kcal = (c.energy - min_energy) * HARTREE_TO_KCAL
                else:
                    # Already in kcal/mol (from MMFF or pre-DFT)
                    energy_kcal = c.energy - min_energy

            conformer_progress.append(
                {
                    "conformer_id": c.conformer_id,
                    "status": c.status,
                    "energy_kcal": round(energy_kcal, 2) if energy_kcal is not None else None,
                    "population": round(c.weight, 4) if c.weight is not None else None,
                }
            )

        # Build ensemble metadata for completed jobs
        if job_status.status == "complete":
            nmr_complete = [c for c in ensemble.conformers if c.status == "nmr_complete"]
            if nmr_complete:
                # Calculate energy range
                energies_kcal = []
                for c in nmr_complete:
                    if c.energy is not None and min_energy is not None:
                        if c.energy_unit == "hartree":
                            energies_kcal.append((c.energy - min_energy) * HARTREE_TO_KCAL)
                        else:
                            energies_kcal.append(c.energy - min_energy)

                energy_range = max(energies_kcal) - min(energies_kcal) if energies_kcal else 0.0

                # Get top 3 by population
                sorted_by_pop = sorted(
                    [c for c in nmr_complete if c.weight is not None],
                    key=lambda x: x.weight,
                    reverse=True,
                )[:3]

                top_populations = []
                for c in sorted_by_pop:
                    e_kcal = None
                    if c.energy is not None and min_energy is not None:
                        if c.energy_unit == "hartree":
                            e_kcal = (c.energy - min_energy) * HARTREE_TO_KCAL
                        else:
                            e_kcal = c.energy - min_energy
                    top_populations.append(
                        {
                            "id": c.conformer_id,
                            "population": round(c.weight, 4),
                            "energy_kcal": round(e_kcal, 2) if e_kcal is not None else None,
                        }
                    )

                ensemble_metadata = {
                    "conformer_count": len(nmr_complete),
                    "total_generated": ensemble.total_generated,
                    "method": ensemble.method,
                    "temperature_k": ensemble.temperature_k,
                    "energy_range_kcal": round(energy_range, 2),
                    "top_populations": top_populations,
                    "conformer_method_warning": conformer_method_warning,
                }

        # Add ensemble_metadata to nmr_results if present
        if nmr_results is not None and ensemble_metadata is not None:
            nmr_results["ensemble_metadata"] = ensemble_metadata

    return {
        "job_id": job_status.job_id,
        "status": job_status.status,
        "created_at": job_status.created_at.isoformat() + "Z",
        "started_at": (
            job_status.started_at.isoformat() + "Z" if job_status.started_at else None
        ),
        "completed_at": (
            job_status.completed_at.isoformat() + "Z"
            if job_status.completed_at
            else None
        ),
        "input_smiles": job_status.input.smiles,
        "input_name": job_status.input.name,
        "preset": job_status.input.preset,
        "solvent": job_status.input.solvent,
        "current_step": job_status.current_step,
        "current_step_label": job_status.current_step_label,
        "step_started_at": (
            job_status.step_started_at.isoformat() + "Z"
            if job_status.step_started_at
            else None
        ),
        "steps_completed": steps_completed,
        "error_message": job_status.error_message,
        "nmr_results": nmr_results,
        "conformer_mode": job_status.conformer_mode,
        # Ensemble-specific fields
        "conformer_method": conformer_method,
        "conformer_count": conformer_count,
        "conformer_progress": conformer_progress,
        "conformer_method_warning": conformer_method_warning,
        "ensemble_metadata": ensemble_metadata,
        # TODO: eta_seconds - implement in Plan 03
        "eta_seconds": None,
    }


@router.post(
    "",
    response_model=JobStatusResponse,
    status_code=status.HTTP_202_ACCEPTED,
    responses={
        202: {"description": "Job accepted for processing"},
        422: {"model": ProblemDetail, "description": "Invalid molecule"},
    },
)
async def submit_smiles(request: JobSubmitRequest):
    """Submit molecule via SMILES string for NMR calculation.

    Returns immediately with job ID. Use GET /api/v1/jobs/{job_id} to poll status.
    """
    # Validate SMILES
    mol, error = validate_smiles(request.smiles)
    if mol is None:
        raise HTTPException(
            status_code=422,
            detail={
                "type": "https://qm-nmr-calc.example/problems/invalid-smiles",
                "title": "Invalid SMILES String",
                "status": 422,
                "detail": error,
            },
        )

    # Validate solvent
    normalized_solvent = validate_solvent(request.solvent)
    if normalized_solvent is None:
        supported = ", ".join(get_supported_solvents())
        raise HTTPException(
            status_code=422,
            detail={
                "type": "https://qm-nmr-calc.example/problems/invalid-solvent",
                "title": "Invalid Solvent",
                "status": 422,
                "detail": f"Unknown solvent '{request.solvent}'. Supported solvents: {supported}",
            },
        )

    # Get software versions for reproducibility
    nwchem_version = get_nwchem_version()

    # Create job directory and initial status
    job_status = create_job_directory(
        smiles=request.smiles,
        solvent=normalized_solvent,
        nwchem_version=nwchem_version,
        name=request.name,
        preset=request.preset,
        notification_email=request.notification_email,
        conformer_mode=request.conformer_mode,
        conformer_method=request.conformer_method,
        max_conformers=request.max_conformers,
    )

    # Generate initial 3D geometry for immediate visualization
    job_dir = get_job_dir(job_status.job_id)
    initial_xyz_path = job_dir / "output" / "initial.xyz"
    _generate_initial_xyz(request.smiles, initial_xyz_path)

    # Queue the appropriate NMR calculation task based on conformer mode
    if job_status.input.conformer_mode == "ensemble":
        run_ensemble_nmr_task(job_status.job_id)
    else:
        run_nmr_task(job_status.job_id)

    return JSONResponse(
        status_code=status.HTTP_202_ACCEPTED,
        content=job_status_to_response(job_status),
        headers={
            "Location": f"/api/v1/jobs/{job_status.job_id}",
            "Retry-After": "30",
        },
    )


@router.post(
    "/upload",
    response_model=JobStatusResponse,
    status_code=status.HTTP_202_ACCEPTED,
    responses={
        202: {"description": "Job accepted for processing"},
        422: {"model": ProblemDetail, "description": "Invalid molecule file"},
    },
)
async def submit_file(
    file: UploadFile,
    solvent: Annotated[
        str, Form(description="NMR solvent for COSMO solvation model")
    ],
    name: Annotated[
        Optional[str], Form(description="Optional molecule name")
    ] = None,
    preset: Annotated[
        str, Form(description="Calculation preset: draft or production")
    ] = "production",
    notification_email: Annotated[
        Optional[str], Form(description="Email for completion notification")
    ] = None,
    conformer_mode: Annotated[
        str, Form(description="Conformational sampling mode: single or ensemble")
    ] = "single",
    conformer_method: Annotated[
        Optional[str], Form(description="Conformer generation method: rdkit_kdg or crest")
    ] = None,
    max_conformers: Annotated[
        Optional[int], Form(description="Maximum conformers to generate")
    ] = None,
):
    """Submit molecule via MOL/SDF file upload for NMR calculation.

    Accepts .mol or .sdf files with a single molecule.
    Returns immediately with job ID. Use GET /api/v1/jobs/{job_id} to poll status.
    """
    # Validate file type by extension
    filename = file.filename or ""
    if not (filename.endswith(".mol") or filename.endswith(".sdf")):
        # Also check by content-type
        allowed_types = {
            "chemical/x-mdl-molfile",
            "chemical/x-mdl-sdfile",
            "application/octet-stream",
        }
        if file.content_type not in allowed_types:
            raise HTTPException(
                status_code=422,
                detail={
                    "type": "https://qm-nmr-calc.example/problems/invalid-file-type",
                    "title": "Invalid File Type",
                    "status": 422,
                    "detail": f"Expected .mol or .sdf file, got '{filename}'",
                },
            )

    # Read and validate content
    content = await file.read()
    mol, error = validate_mol_file(content, filename)
    if mol is None:
        raise HTTPException(
            status_code=422,
            detail={
                "type": "https://qm-nmr-calc.example/problems/invalid-molecule-file",
                "title": "Invalid Molecule File",
                "status": 422,
                "detail": error,
            },
        )

    # Validate solvent
    normalized_solvent = validate_solvent(solvent)
    if normalized_solvent is None:
        supported = ", ".join(get_supported_solvents())
        raise HTTPException(
            status_code=422,
            detail={
                "type": "https://qm-nmr-calc.example/problems/invalid-solvent",
                "title": "Invalid Solvent",
                "status": 422,
                "detail": f"Unknown solvent '{solvent}'. Supported solvents: {supported}",
            },
        )

    # Convert to SMILES for storage
    smiles = Chem.MolToSmiles(mol)

    # Get software versions
    nwchem_version = get_nwchem_version()

    # Create job (use filename as name if not provided)
    job_status = create_job_directory(
        smiles=smiles,
        solvent=normalized_solvent,
        nwchem_version=nwchem_version,
        name=name or filename,
        preset=preset,
        notification_email=notification_email,
        conformer_mode=conformer_mode,
        conformer_method=conformer_method,
        max_conformers=max_conformers,
    )

    # Generate initial 3D geometry for immediate visualization
    job_dir = get_job_dir(job_status.job_id)
    initial_xyz_path = job_dir / "output" / "initial.xyz"
    _generate_initial_xyz(smiles, initial_xyz_path)

    # Queue the appropriate NMR calculation task based on conformer mode
    if job_status.input.conformer_mode == "ensemble":
        run_ensemble_nmr_task(job_status.job_id)
    else:
        run_nmr_task(job_status.job_id)

    return JSONResponse(
        status_code=status.HTTP_202_ACCEPTED,
        content=job_status_to_response(job_status),
        headers={
            "Location": f"/api/v1/jobs/{job_status.job_id}",
            "Retry-After": "30",
        },
    )


@router.get(
    "/solvents",
    response_model=list[str],
    responses={
        200: {"description": "List of supported solvents"},
    },
)
async def list_solvents():
    """List supported NMR solvents for COSMO solvation model.

    Returns list of NWChem COSMO solvent names that can be used
    in job submission requests.
    """
    return get_supported_solvents()


@router.get(
    "/{job_id}",
    response_model=JobStatusResponse,
    responses={
        200: {"description": "Job status"},
        404: {"model": ProblemDetail, "description": "Job not found"},
    },
)
async def get_job_status(job_id: str):
    """Get status of a calculation job.

    Returns full job status including timestamps and any error message.
    """
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )
    return job_status_to_response(job_status)


@router.get(
    "/{job_id}/results",
    response_model=NMRResultsResponse,
    responses={
        200: {"description": "NMR calculation results"},
        404: {"model": ProblemDetail, "description": "Job or results not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def get_nmr_results(job_id: str):
    """Get NMR chemical shifts for a completed job.

    Returns 1H and 13C shifts with atom indices, plus calculation metadata.
    """
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )

    if job_status.status != "complete":
        raise HTTPException(
            status_code=409,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-complete",
                "title": "Job Not Complete",
                "status": 409,
                "detail": f"Job '{job_id}' is in '{job_status.status}' state. Results available when complete.",
            },
        )

    if job_status.nmr_results is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/results-not-found",
                "title": "Results Not Found",
                "status": 404,
                "detail": f"Job '{job_id}' completed but no NMR results available.",
            },
        )

    # Get factor metadata for MAE display
    # Note: functional stored lowercase, get_scaling_factor expects uppercase
    h1_factor = get_scaling_factor(
        job_status.nmr_results.functional.upper(),
        job_status.nmr_results.basis_set,
        "1H",
        job_status.nmr_results.solvent,
    )
    c13_factor = get_scaling_factor(
        job_status.nmr_results.functional.upper(),
        job_status.nmr_results.basis_set,
        "13C",
        job_status.nmr_results.solvent,
    )

    # Build ensemble metadata if this is an ensemble job
    ensemble_metadata = None
    if job_status.conformer_mode == "ensemble" and job_status.conformer_ensemble:
        from ...conformers.boltzmann import HARTREE_TO_KCAL

        ensemble = job_status.conformer_ensemble
        conformer_method_warning = getattr(job_status, "conformer_method_warning", None)

        # Find minimum energy for relative energy display
        min_energy = None
        for c in ensemble.conformers:
            if c.energy is not None and c.status == "nmr_complete":
                if min_energy is None or c.energy < min_energy:
                    min_energy = c.energy

        nmr_complete = [c for c in ensemble.conformers if c.status == "nmr_complete"]
        if nmr_complete:
            # Calculate energy range
            energies_kcal = []
            for c in nmr_complete:
                if c.energy is not None and min_energy is not None:
                    if c.energy_unit == "hartree":
                        energies_kcal.append((c.energy - min_energy) * HARTREE_TO_KCAL)
                    else:
                        energies_kcal.append(c.energy - min_energy)

            energy_range = max(energies_kcal) - min(energies_kcal) if energies_kcal else 0.0

            # Get top 3 by population
            sorted_by_pop = sorted(
                [c for c in nmr_complete if c.weight is not None],
                key=lambda x: x.weight,
                reverse=True,
            )[:3]

            top_populations = []
            for c in sorted_by_pop:
                e_kcal = None
                if c.energy is not None and min_energy is not None:
                    if c.energy_unit == "hartree":
                        e_kcal = (c.energy - min_energy) * HARTREE_TO_KCAL
                    else:
                        e_kcal = c.energy - min_energy
                top_populations.append(
                    {
                        "id": c.conformer_id,
                        "population": round(c.weight, 4),
                        "energy_kcal": round(e_kcal, 2) if e_kcal is not None else None,
                    }
                )

            from ..schemas import EnsembleMetadataResponse

            ensemble_metadata = EnsembleMetadataResponse(
                conformer_count=len(nmr_complete),
                total_generated=ensemble.total_generated,
                method=ensemble.method,
                temperature_k=ensemble.temperature_k,
                energy_range_kcal=round(energy_range, 2),
                top_populations=top_populations,
                conformer_method_warning=conformer_method_warning,
            )

    # Convert to response format (without shielding, only shift)
    return NMRResultsResponse(
        h1_shifts=[
            {"index": s.index, "atom": s.atom, "shift": s.shift}
            for s in job_status.nmr_results.h1_shifts
        ],
        c13_shifts=[
            {"index": s.index, "atom": s.atom, "shift": s.shift}
            for s in job_status.nmr_results.c13_shifts
        ],
        functional=job_status.nmr_results.functional,
        basis_set=job_status.nmr_results.basis_set,
        solvent=job_status.nmr_results.solvent,
        scaling_factor_source="DELTA50",
        h1_expected_mae=f"+/- {h1_factor['mae']:.2f} ppm",
        c13_expected_mae=f"+/- {c13_factor['mae']:.2f} ppm",
        ensemble_metadata=ensemble_metadata,
    )


@router.get(
    "/{job_id}/geometry",
    response_class=FileResponse,
    responses={
        200: {"description": "Optimized molecular geometry (XYZ)", "content": {"chemical/x-xyz": {}}},
        404: {"model": ProblemDetail, "description": "Job or geometry file not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_geometry(job_id: str):
    """Download optimized molecular geometry as XYZ file."""
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )

    if job_status.status != "complete":
        raise HTTPException(
            status_code=409,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-complete",
                "title": "Job Not Complete",
                "status": 409,
                "detail": f"Job '{job_id}' is in '{job_status.status}' state. Geometry available when complete.",
            },
        )

    geometry_file = get_geometry_file(job_id)
    if geometry_file is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/geometry-not-found",
                "title": "Geometry File Not Found",
                "status": 404,
                "detail": f"Optimized geometry file not found for job '{job_id}'.",
            },
        )

    return FileResponse(
        path=geometry_file,
        media_type="chemical/x-xyz",
        filename=f"{job_id}_optimized.xyz",
    )


@router.get(
    "/{job_id}/geometry.sdf",
    response_class=Response,
    responses={
        200: {"description": "Optimized molecular geometry (SDF)", "content": {"chemical/x-mdl-sdfile": {}}},
        404: {"model": ProblemDetail, "description": "Job or geometry file not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_geometry_sdf(job_id: str):
    """Download optimized molecular geometry as SDF file.

    Generates SDF from original SMILES with optimized 3D coordinates.
    """
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )

    if job_status.status != "complete":
        raise HTTPException(
            status_code=409,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-complete",
                "title": "Job Not Complete",
                "status": 409,
                "detail": f"Job '{job_id}' is in '{job_status.status}' state. Geometry available when complete.",
            },
        )

    geometry_file = get_geometry_file(job_id)
    if geometry_file is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/geometry-not-found",
                "title": "Geometry File Not Found",
                "status": 404,
                "detail": f"Optimized geometry file not found for job '{job_id}'.",
            },
        )

    # Read XYZ coordinates
    xyz_content = geometry_file.read_text()
    xyz_lines = xyz_content.strip().split("\n")

    # Parse XYZ: first line is atom count, second is comment, rest are coords
    coords = []
    for line in xyz_lines[2:]:  # Skip count and comment lines
        parts = line.split()
        if len(parts) >= 4:
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

    # Create mol from SMILES and set coordinates
    mol = Chem.MolFromSmiles(job_status.input.smiles)
    if mol is None:
        raise HTTPException(
            status_code=500,
            detail={
                "type": "https://qm-nmr-calc.example/problems/sdf-generation-failed",
                "title": "SDF Generation Failed",
                "status": 500,
                "detail": "Failed to parse original SMILES for SDF generation.",
            },
        )

    mol = Chem.AddHs(mol)

    # Create conformer with coordinates
    from rdkit.Chem import AllChem
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        if i < mol.GetNumAtoms():
            conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf, assignId=True)

    # Generate SDF content
    sdf_content = Chem.MolToMolBlock(mol)

    return Response(
        content=sdf_content,
        media_type="chemical/x-mdl-sdfile",
        headers={"Content-Disposition": f'attachment; filename="{job_id}_optimized.sdf"'},
    )


@router.get(
    "/{job_id}/output",
    responses={
        200: {"description": "Raw NWChem output files (ZIP)", "content": {"application/zip": {}}},
        404: {"model": ProblemDetail, "description": "Job or output files not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_output(job_id: str):
    """Download raw NWChem output files as ZIP archive."""
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )

    if job_status.status != "complete":
        raise HTTPException(
            status_code=409,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-complete",
                "title": "Job Not Complete",
                "status": 409,
                "detail": f"Job '{job_id}' is in '{job_status.status}' state. Output files available when complete.",
            },
        )

    output_files = get_output_files(job_id)
    if not output_files:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/output-not-found",
                "title": "Output Files Not Found",
                "status": 404,
                "detail": f"No NWChem output files found for job '{job_id}'.",
            },
        )

    # Create in-memory ZIP
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        for filepath in output_files:
            zf.write(filepath, filepath.name)
    zip_buffer.seek(0)

    return Response(
        content=zip_buffer.getvalue(),
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="{job_id}_output.zip"'},
    )


async def _get_visualization(job_id: str, filename: str, media_type: str, download_name: str):
    """Common logic for visualization download endpoints."""
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )

    if job_status.status != "complete":
        raise HTTPException(
            status_code=409,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-complete",
                "title": "Job Not Complete",
                "status": 409,
                "detail": f"Job '{job_id}' is in '{job_status.status}' state. Visualizations available when complete.",
            },
        )

    viz_file = get_visualization_file(job_id, filename)
    if viz_file is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/visualization-not-found",
                "title": "Visualization Not Found",
                "status": 404,
                "detail": f"Visualization file '{filename}' not found for job '{job_id}'.",
            },
        )

    return FileResponse(
        path=viz_file,
        media_type=media_type,
        filename=download_name,
    )


@router.get(
    "/{job_id}/spectrum/1h.png",
    response_class=FileResponse,
    responses={
        200: {"description": "1H NMR spectrum plot (PNG)", "content": {"image/png": {}}},
        404: {"model": ProblemDetail, "description": "Job or file not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_h1_spectrum_png(job_id: str):
    """Download 1H NMR spectrum plot as PNG (300 DPI)."""
    return await _get_visualization(job_id, "spectrum_1H.png", "image/png", f"{job_id}_1H_spectrum.png")


@router.get(
    "/{job_id}/spectrum/1h.svg",
    response_class=FileResponse,
    responses={
        200: {"description": "1H NMR spectrum plot (SVG)", "content": {"image/svg+xml": {}}},
        404: {"model": ProblemDetail, "description": "Job or file not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_h1_spectrum_svg(job_id: str):
    """Download 1H NMR spectrum plot as SVG."""
    return await _get_visualization(job_id, "spectrum_1H.svg", "image/svg+xml", f"{job_id}_1H_spectrum.svg")


@router.get(
    "/{job_id}/spectrum/13c.png",
    response_class=FileResponse,
    responses={
        200: {"description": "13C NMR spectrum plot (PNG)", "content": {"image/png": {}}},
        404: {"model": ProblemDetail, "description": "Job or file not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_c13_spectrum_png(job_id: str):
    """Download 13C NMR spectrum plot as PNG (300 DPI)."""
    return await _get_visualization(job_id, "spectrum_13C.png", "image/png", f"{job_id}_13C_spectrum.png")


@router.get(
    "/{job_id}/spectrum/13c.svg",
    response_class=FileResponse,
    responses={
        200: {"description": "13C NMR spectrum plot (SVG)", "content": {"image/svg+xml": {}}},
        404: {"model": ProblemDetail, "description": "Job or file not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_c13_spectrum_svg(job_id: str):
    """Download 13C NMR spectrum plot as SVG."""
    return await _get_visualization(job_id, "spectrum_13C.svg", "image/svg+xml", f"{job_id}_13C_spectrum.svg")


@router.get(
    "/{job_id}/structure.png",
    response_class=FileResponse,
    responses={
        200: {"description": "Annotated structure image (PNG)", "content": {"image/png": {}}},
        404: {"model": ProblemDetail, "description": "Job or file not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_structure_png(job_id: str):
    """Download annotated molecular structure as PNG (300 DPI equivalent)."""
    return await _get_visualization(job_id, "structure_annotated.png", "image/png", f"{job_id}_structure.png")


@router.get(
    "/{job_id}/structure.svg",
    response_class=FileResponse,
    responses={
        200: {"description": "Annotated structure image (SVG)", "content": {"image/svg+xml": {}}},
        404: {"model": ProblemDetail, "description": "Job or file not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_structure_svg(job_id: str):
    """Download annotated molecular structure as SVG."""
    return await _get_visualization(job_id, "structure_annotated.svg", "image/svg+xml", f"{job_id}_structure.svg")


def _xyz_to_sdf(xyz_content: str, smiles: str) -> str | None:
    """Convert XYZ content to SDF using SMILES for bond information.

    Args:
        xyz_content: XYZ format geometry string
        smiles: SMILES string for bond connectivity

    Returns:
        SDF/MOL block string, or None if conversion fails
    """
    try:
        xyz_lines = xyz_content.strip().split("\n")
        coords = []
        for line in xyz_lines[2:]:  # Skip count and comment lines
            parts = line.split()
            if len(parts) >= 4:
                coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = Chem.AddHs(mol)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            if i < mol.GetNumAtoms():
                conf.SetAtomPosition(i, (x, y, z))
        mol.AddConformer(conf, assignId=True)
        return Chem.MolToMolBlock(mol)
    except Exception:
        return None


@router.get(
    "/{job_id}/geometry.json",
    responses={
        200: {"description": "Geometry and shift data for 3D visualization"},
        404: {"model": ProblemDetail, "description": "Job not found"},
    },
)
async def get_geometry_data(job_id: str):
    """Get geometry and shift data for 3D visualization.

    For single-conformer jobs:
        Returns initial RDKit geometry for running jobs,
        optimized NWChem geometry for complete jobs.

    For ensemble jobs:
        Returns array of conformer geometries with energies and populations.
        Default geometry (xyz/sdf) is the lowest-energy conformer.
    """
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )

    # Build shift assignments only for complete jobs
    h1_assignments = {}
    c13_assignments = {}
    if job_status.status == "complete" and job_status.nmr_results:
        for s in job_status.nmr_results.h1_shifts:
            h1_assignments[str(s.index)] = s.shift
        for s in job_status.nmr_results.c13_shifts:
            c13_assignments[str(s.index)] = s.shift

    # Handle ensemble mode
    if job_status.conformer_mode == "ensemble" and job_status.conformer_ensemble:
        from ...conformers.boltzmann import HARTREE_TO_KCAL

        ensemble = job_status.conformer_ensemble
        job_dir = get_job_dir(job_id)

        # Build conformer array
        conformers_data = []
        min_energy = None

        # Find minimum energy for relative calculation
        for c in ensemble.conformers:
            if c.energy is not None and c.status in ("optimized", "nmr_complete"):
                if min_energy is None or c.energy < min_energy:
                    min_energy = c.energy

        # Sort by energy (lowest first)
        sorted_conformers = sorted(
            [c for c in ensemble.conformers if c.status in ("optimized", "nmr_complete")],
            key=lambda x: x.energy if x.energy is not None else float("inf"),
        )

        for c in sorted_conformers:
            # Read geometry file (prefer optimized over initial)
            geom_file = c.optimized_geometry_file or c.geometry_file
            if geom_file:
                geom_path = job_dir / geom_file
                if geom_path.exists():
                    xyz_content = geom_path.read_text()

                    # Calculate relative energy in kcal/mol
                    energy_kcal = None
                    if c.energy is not None and min_energy is not None:
                        if c.energy_unit == "hartree":
                            energy_kcal = (c.energy - min_energy) * HARTREE_TO_KCAL
                        else:
                            energy_kcal = c.energy - min_energy

                    # Generate SDF from XYZ + SMILES for proper bonds
                    sdf_content = None
                    if job_status.status == "complete" and job_status.input.smiles:
                        sdf_content = _xyz_to_sdf(xyz_content, job_status.input.smiles)

                    conformers_data.append({
                        "id": c.conformer_id,
                        "xyz": xyz_content,
                        "sdf": sdf_content,
                        "energy_kcal": round(energy_kcal, 2) if energy_kcal is not None else None,
                        "population": round(c.weight, 4) if c.weight is not None else None,
                    })

        # Default geometry is first conformer (lowest energy)
        default_xyz = conformers_data[0]["xyz"] if conformers_data else None
        default_sdf = conformers_data[0]["sdf"] if conformers_data else None

        return {
            "job_id": job_id,
            "status": job_status.status,
            "conformer_mode": "ensemble",
            "xyz": default_xyz,
            "sdf": default_sdf,
            "h1_assignments": h1_assignments if h1_assignments else None,
            "c13_assignments": c13_assignments if c13_assignments else None,
            "conformers": conformers_data,
        }

    # Single conformer mode (existing logic)
    # Determine which geometry to return
    # Prefer optimized geometry if available (even during NMR shielding step)
    geometry_file = get_geometry_file(job_id)  # optimized.xyz
    if geometry_file is None:
        geometry_file = get_initial_geometry_file(job_id)  # initial.xyz fallback

    if geometry_file is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/geometry-not-found",
                "title": "Geometry Not Found",
                "status": 404,
                "detail": f"Geometry file not available for job '{job_id}'",
            },
        )

    xyz_content = geometry_file.read_text()

    # Generate SDF for complete jobs
    sdf_content = None
    if job_status.status == "complete" and job_status.nmr_results and job_status.input.smiles:
        sdf_content = _xyz_to_sdf(xyz_content, job_status.input.smiles)

    return {
        "job_id": job_id,
        "status": job_status.status,
        "conformer_mode": "single",
        "xyz": xyz_content,
        "sdf": sdf_content,
        "h1_assignments": h1_assignments if h1_assignments else None,
        "c13_assignments": c13_assignments if c13_assignments else None,
        "conformers": None,
    }
